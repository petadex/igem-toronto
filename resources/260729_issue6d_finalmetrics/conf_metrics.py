"""
conf_metrics.py -- figure Da: the self-confidence distribution of a finished fold run.

Reads mean pLDDT / pTM for every folded ORF in one prediction run and answers "how many of these
structures are actually usable?" -- no ground truth, no second run, no structure parsing.

    python conf_metrics.py --run s3://bucket/esmfold2-centroids/ --out conf_centroids

`--run` is an `s3://bucket/prefix/` or a local directory laid out the way
esmfold2_local_predictor.py writes one:  <root>/metrics/orf<ORFid>.json.

WHY metrics/*.json AND NOT tracker.csv
  tracker.csv carries the same scalars in one small file, which looks like the cheap way in. It is
  not the trustworthy one at scale:
    * the full centroid run is sharded across a fleet, and each worker mirrors its OWN
      tracker.shardKofN.csv -- so "the ledger" is really N files that have to be found and unioned;
    * each mirror happens on a ~300 s interval (MIRROR_EVERY_S), so a tracker read mid-run is stale
      by construction, and a replaced spot instance can come back with a ledger BEHIND the objects
      it already uploaded.
  A metrics/<id>.json object, by contrast, exists if and only if that fold finished AND uploaded.
  The listing IS the set of structures we have, which is exactly the population figure Da is about.

WHY IT IS STILL CHEAP (this is the part that makes the choice free)
  The rich per-protein record is dominated by `pae`, which is L x L -- a 600 aa ORF's metrics JSON is
  ~2 MB of text, so GETting them whole at 630k would move ~600 GB. But the scalars we need
  (`seq_len`, `confidence.mean_plddt`, `confidence.ptm`) are written BEFORE `per_residue_plddt` and
  `pae` in the object, so a Range-GET of the first few KB has all of them:

      630k x 4 KB ~= 2.5 GB transferred, 630k GET requests ~= $0.25

  The run is request-latency-bound, not bandwidth-bound, so --workers scales it nearly linearly.
  --full-get forces whole objects if the layout ever changes; --range-bytes widens the window.

BUILT FOR THE FULL CENTROID SET, same conventions as pair_tm.py:
  * `--shard K/N` fans one command across a fleet; each shard owns a disjoint id band and its own
    ledger, so workers never clobber each other;
  * the ledger is append-only and re-read on start, so an interrupted run RESUMES;
  * per-ORF rows are a convenience -- the plotting input is `aggregates*.json` (fixed-bin histograms
    + Welford moments + EXACT threshold counts + a reservoir subsample), which is O(1) in N and
    mergeable across shards (`--merge`). The figures cost the same at 121 ORFs or 2 million.

READ THIS BEFORE INTERPRETING:
  * The denominator here is STRUCTURES THAT EXIST, not ORFs attempted. Folds that errored never
    wrote a metrics JSON, so they are invisible to this script -- by design, since Da asks how usable
    the structures we HAVE are. The attempt-level failure rate is a different question and must come
    from the input id list, not from here.
  * pLDDT/pTM are the model's own confidence, i.e. a PROXY for quality. High confidence is not
    proof of correctness; it is the standard stand-in when ground truth is unavailable at scale
    (which is the whole premise of group D). Ground-truth agreement is group C.
  * `truncated` rows are ORFs longer than the ESMC context (2048 aa) whose tail was cut before
    folding. Their confidence describes the folded PREFIX, not the ORF, so they are counted and
    reported separately rather than silently pooled.

Deps:  pip install numpy boto3
"""
import argparse
import csv
import json
import math
import os
import random
import re
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import numpy as np


# ---------------------------------------------------------------------------
# reading the scalars out of the head of a metrics JSON
# ---------------------------------------------------------------------------
# The head is a TRUNCATED JSON document, so it cannot be json.loads()'d -- these pull the scalars out
# of the raw text instead. Each pattern requires the opening quote of the key, so `"ptm"` does NOT
# match inside `"iptm"` (the char before `ptm` there is `i`, not `"`). Keys are unique in the head and
# array contents are bare numbers, never keys, so there is nothing else for these to hit.
_NUM = r"(-?[\d.eE+]+)"
_RX = {
    "seq_len":     re.compile(r'"seq_len":\s*(\d+)'),
    "truncated":   re.compile(r'"truncated":\s*(true|false)'),
    "wall_s":      re.compile(r'"wall_s_median":\s*' + _NUM),
    "res_per_s":   re.compile(r'"residues_per_s":\s*' + _NUM),
    "mean_plddt":  re.compile(r'"mean_plddt":\s*' + _NUM),
    "ptm":         re.compile(r'"ptm":\s*' + _NUM),
    "iptm":        re.compile(r'"iptm":\s*' + _NUM),
}
# Present in every record, so their absence means the head window was too small (or the layout
# changed) -- the fetcher retries such an object with a full GET rather than recording a bad row.
_REQUIRED = ("mean_plddt", "ptm")


def parse_head(text):
    """Scalars from the (possibly truncated) head of a metrics JSON. Missing keys -> None."""
    out = {}
    for k, rx in _RX.items():
        m = rx.search(text)
        if not m:
            out[k] = None
            continue
        v = m.group(1)
        if k == "truncated":
            out[k] = v == "true"
        elif k == "seq_len":
            out[k] = int(v)
        else:
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = None
    return out


class MetricSource:
    """Reads <root>/metrics/<json_name> as text. `root` is s3://bucket/prefix/ or a local dir.

    One boto3 client PER THREAD: clients are only nominally thread-safe, and sharing one across a
    large pool serialises on its connection pool (the whole point here is concurrent GETs).
    """

    def __init__(self, root, json_name="orf{id}.json", range_bytes=4096, full_get=False, retries=4):
        self.root, self.json_name, self.retries = root, json_name, retries
        self.range_bytes, self.full_get = range_bytes, full_get
        self.is_s3 = str(root).startswith("s3://")
        self._tl = threading.local()
        self.n_widened = 0                                # objects that needed a full GET (diagnostic)
        if self.is_s3:
            self.bucket, _, self.prefix = root[len("s3://"):].partition("/")
            self.prefix = self.prefix.strip("/")
        self.label = os.path.basename(str(root).rstrip("/\\")) or str(root)

    @property
    def s3(self):
        if getattr(self._tl, "c", None) is None:
            import boto3                                  # pre-warmed by warm_imports()
            self._tl.c = boto3.client("s3")
        return self._tl.c

    def _key(self, oid):
        return "/".join(p for p in (self.prefix, "metrics", self.json_name.format(id=oid)) if p)

    def _get(self, oid, whole):
        """Object text (head slice unless `whole`), or None if it is not there."""
        name = self.json_name.format(id=oid)
        if not self.is_s3:
            for p in (os.path.join(self.root, "metrics", name), os.path.join(self.root, name)):
                if os.path.exists(p):
                    with open(p) as f:
                        return f.read() if whole else f.read(self.range_bytes)
            return None
        kw = {"Bucket": self.bucket, "Key": self._key(oid)}
        if not whole:
            kw["Range"] = f"bytes=0-{self.range_bytes - 1}"
        for attempt in range(self.retries):
            try:
                return self.s3.get_object(**kw)["Body"].read().decode()
            except Exception as e:                        # noqa: BLE001
                code = str(getattr(e, "response", {}).get("Error", {}).get("Code", ""))
                if code in ("NoSuchKey", "404", "NoSuchBucket") or type(e).__name__ == "NoSuchKey":
                    return None
                if attempt == self.retries - 1:
                    raise
                time.sleep(0.4 * 2 ** attempt)            # throttling / transient -> backoff
        return None

    def scalars(self, oid):
        """(record dict, status). Falls back to a full GET if the head lacked a required scalar."""
        text = self._get(oid, whole=self.full_get)
        if text is None:
            return {}, "missing"
        rec = parse_head(text)
        if not self.full_get and any(rec.get(k) is None for k in _REQUIRED):
            text = self._get(oid, whole=True)             # window too small / layout changed
            if text is None:
                return {}, "missing"
            self.n_widened += 1
            rec = parse_head(text)
        if any(rec.get(k) is None for k in _REQUIRED):
            return rec, "parse_fail"
        return rec, "ok"

    def list_ids(self):
        """Every ORFid with a metrics JSON under this root (paginated; 1000 keys/request).

        This listing IS the work list: an object exists iff that fold finished and uploaded.
        """
        pat = re.compile("^" + re.escape(self.json_name).replace(r"\{id\}", "(.+)") + "$")
        out = set()
        if not self.is_s3:
            d = os.path.join(self.root, "metrics")
            for n in (os.listdir(d) if os.path.isdir(d) else []):
                m = pat.match(n)
                if m:
                    out.add(norm_id(m.group(1)))
            return out
        tok = None
        while True:
            kw = {"Bucket": self.bucket, "Prefix": "/".join(p for p in (self.prefix, "metrics/") if p)}
            if tok:
                kw["ContinuationToken"] = tok
            r = self.s3.list_objects_v2(**kw)
            for o in r.get("Contents", []):
                m = pat.match(o["Key"].rsplit("/", 1)[-1])
                if m:
                    out.add(norm_id(m.group(1)))
            if not r.get("IsTruncated"):
                return out
            tok = r["NextContinuationToken"]


def norm_id(x):
    """Canonical bare ORFid. Files are named orf<ORFid>.json while map files use bare <ORFid>; both
    must land on the same key so ledgers and joins line up."""
    x = str(x).strip()
    return x[3:] if x.startswith("orf") else x


def warm_imports(needs_s3):
    """Import boto3 ONCE on the main thread. Concurrent first-imports inside a pool race on
    sys.modules and intermittently blow up with half-initialised module errors."""
    if needs_s3:
        import boto3  # noqa: F401


# ---------------------------------------------------------------------------
# sharding (same 'K/N' contract as esmfold2_local_predictor.py --shard)
# ---------------------------------------------------------------------------
def parse_shard(spec):
    if not spec:
        return None
    try:
        ks, ns = str(spec).split("/")
        k, n = int(ks), int(ns)
    except Exception:                                     # noqa: BLE001
        raise SystemExit(f"--shard must be 'K/N' (1-based), got {spec!r}")
    if not (1 <= k <= n):
        raise SystemExit(f"--shard 'K/N' needs 1<=K<=N, got {k}/{n}")
    return k, n


def select_shard(ids, k, n):
    """Contiguous band K of N over deterministically sorted ids -> shards never overlap or gap,
    whatever order the listing returned them in."""
    ordered = sorted(ids, key=lambda s: (len(s), s))
    b = [round(i * len(ordered) / n) for i in range(n + 1)]
    return ordered[b[k - 1]:b[k]]


# ---------------------------------------------------------------------------
# streaming aggregates -- the actual plotting input. O(1) in N, mergeable across shards.
# ---------------------------------------------------------------------------
# Histogram edges are FIXED (never data-derived) so shard JSONs stay bin-for-bin addable.
PLDDT_BINS = (0.0, 100.0, 200)      # 0.5 pLDDT wide
PTM_BINS = (0.0, 1.0, 200)          # 0.005 wide
LEN_BINS = (0.0, 2100.0, 42)        # 50 aa wide (ESMC context is 2048)
JOINT_NB = 100                      # 2D panels: 100 x 100 cells
# Da names 70 and 80; the rest bracket them so any threshold the writeup lands on is already exact.
PLDDT_THRESHOLDS = (50.0, 60.0, 70.0, 80.0, 90.0)
PTM_THRESHOLDS = (0.5, 0.7, 0.8, 0.9)


def _edges(spec):
    lo, hi, nb = spec
    return np.linspace(lo, hi, nb + 1)


def _bin(v, spec):
    """Index into `spec`'s bins, with everything >= hi landing in one extra overflow cell."""
    lo, hi, nb = spec
    if v is None or not (v == v):                         # None / NaN
        return nb
    if v >= hi:
        return nb
    return max(0, min(nb - 1, int((v - lo) / (hi - lo) * nb)))


class _Moments:
    """Welford accumulator: exact streaming mean/variance, and combinable across shards."""

    def __init__(self):
        self.n, self.mean, self.m2 = 0, 0.0, 0.0
        self.min, self.max = float("inf"), float("-inf")

    def add(self, x):
        self.n += 1
        d = x - self.mean
        self.mean += d / self.n
        self.m2 += d * (x - self.mean)
        self.min, self.max = min(self.min, x), max(self.max, x)

    def to_dict(self):
        return {"n": self.n, "mean": self.mean, "m2": self.m2,
                "var": (self.m2 / (self.n - 1)) if self.n > 1 else 0.0,
                "min": self.min if self.n else None, "max": self.max if self.n else None}

    @staticmethod
    def merge_dicts(ds):
        n, mean, m2 = 0, 0.0, 0.0
        mn, mx = float("inf"), float("-inf")
        for d in ds:
            dn = d["n"]
            if not dn:
                continue
            delta = d["mean"] - mean
            tot = n + dn
            mean += delta * dn / tot
            m2 += d["m2"] + delta ** 2 * n * dn / tot
            n = tot
            mn, mx = min(mn, d["min"]), max(mx, d["max"])
        return {"n": n, "mean": mean, "m2": m2, "var": (m2 / (n - 1)) if n > 1 else 0.0,
                "min": mn if n else None, "max": mx if n else None}


class Aggregates:
    """Fixed-bin histograms + Welford moments + EXACT threshold counts + a reservoir subsample.

    Threshold counts are accumulated directly rather than read off the histograms: "fraction above
    70" is the headline number of figure Da, and a bin-edge answer to it would be off by up to a
    bin. The histograms are for shape only.
    """

    def __init__(self, sample_size=20000, seed=0):
        self.n = 0
        self.plddt = np.zeros(PLDDT_BINS[2] + 1, dtype=np.int64)
        self.ptm = np.zeros(PTM_BINS[2] + 1, dtype=np.int64)
        self.length = np.zeros(LEN_BINS[2] + 1, dtype=np.int64)
        self.plddt_by_len = np.zeros((LEN_BINS[2] + 1, JOINT_NB), dtype=np.int64)
        self.plddt_by_ptm = np.zeros((JOINT_NB, JOINT_NB), dtype=np.int64)   # rows pLDDT, cols pTM
        self.plddt_above = {t: 0 for t in PLDDT_THRESHOLDS}
        self.ptm_above = {t: 0 for t in PTM_THRESHOLDS}
        self.m_plddt, self.m_ptm = _Moments(), _Moments()
        self.status = Counter()
        self.n_truncated = 0                              # folded prefix only -- reported separately
        self.trunc_plddt_above = {t: 0 for t in PLDDT_THRESHOLDS}
        self.sample_size = sample_size
        self.sample = []
        self._seen = 0
        self._rng = random.Random(seed)

    # -- one finished ORF --------------------------------------------------
    def add(self, rec):
        self.status[rec["status"]] += 1
        if rec["status"] != "ok":
            return
        pl, pt = float(rec["mean_plddt"]), float(rec["ptm"])
        L = rec.get("seq_len") or 0
        self.n += 1
        self.plddt[_bin(pl, PLDDT_BINS)] += 1
        self.ptm[_bin(pt, PTM_BINS)] += 1
        self.length[_bin(float(L), LEN_BINS)] += 1
        self.plddt_by_len[_bin(float(L), LEN_BINS), min(int(pl), JOINT_NB - 1)] += 1
        self.plddt_by_ptm[min(int(pl), JOINT_NB - 1), min(int(pt * JOINT_NB), JOINT_NB - 1)] += 1
        for t in PLDDT_THRESHOLDS:
            if pl >= t:
                self.plddt_above[t] += 1
        for t in PTM_THRESHOLDS:
            if pt >= t:
                self.ptm_above[t] += 1
        self.m_plddt.add(pl)
        self.m_ptm.add(pt)
        if rec.get("truncated") in (True, "True", "true", 1, "1"):
            self.n_truncated += 1
            for t in PLDDT_THRESHOLDS:
                if pl >= t:
                    self.trunc_plddt_above[t] += 1
        self._reservoir(rec)

    def _reservoir(self, rec):
        """Algorithm R, seeded -> a uniform sample of the stream for the scatter panels at any N,
        reproducible across reruns of the same shard."""
        self._seen += 1
        row = {"orfid": rec["orfid"], "mean_plddt": round(float(rec["mean_plddt"]), 3),
               "ptm": round(float(rec["ptm"]), 4), "seq_len": rec.get("seq_len") or 0}
        if len(self.sample) < self.sample_size:
            self.sample.append(row)
        else:
            j = self._rng.randrange(self._seen)
            if j < self.sample_size:
                self.sample[j] = row

    # -- serialise / merge -------------------------------------------------
    def to_dict(self):
        return {
            "n_ok": self.n, "status": dict(self.status), "n_truncated": self.n_truncated,
            "bins": {"plddt": PLDDT_BINS, "ptm": PTM_BINS, "len": LEN_BINS, "joint_nb": JOINT_NB},
            "hist": {"plddt": self.plddt.tolist(), "ptm": self.ptm.tolist(),
                     "len": self.length.tolist(), "plddt_by_len": self.plddt_by_len.tolist(),
                     "plddt_by_ptm": self.plddt_by_ptm.tolist()},
            "plddt_above": {str(k): v for k, v in self.plddt_above.items()},
            "ptm_above": {str(k): v for k, v in self.ptm_above.items()},
            "trunc_plddt_above": {str(k): v for k, v in self.trunc_plddt_above.items()},
            "moments": {"plddt": self.m_plddt.to_dict(), "ptm": self.m_ptm.to_dict()},
            "sample": {"size": self.sample_size, "seen": self._seen, "rows": self.sample},
        }

    @staticmethod
    def merge(dicts, seed=0):
        """Combine shard aggregates into one. Histograms and exact counts add; Welford combines by
        the parallel formula; reservoirs are resampled from the union weighted by each shard's stream
        length (each shard's sample is uniform over its own stream, so this is uniform over the
        union)."""
        dicts = [d for d in dicts if d]
        if not dicts:
            raise SystemExit("nothing to merge")
        out = {"n_ok": 0, "status": Counter(), "n_truncated": 0, "bins": dicts[0]["bins"],
               "hist": {k: np.zeros_like(np.asarray(v, dtype=np.int64))
                        for k, v in dicts[0]["hist"].items()}}
        counters = {k: Counter() for k in ("plddt_above", "ptm_above", "trunc_plddt_above")}
        for d in dicts:
            if d["bins"] != out["bins"]:
                raise SystemExit("bin edges differ between shards -- cannot merge")
            out["n_ok"] += d["n_ok"]
            out["n_truncated"] += d.get("n_truncated", 0)
            out["status"].update(d["status"])
            for k in out["hist"]:
                out["hist"][k] += np.asarray(d["hist"][k], dtype=np.int64)
            for k, c in counters.items():
                c.update(d.get(k, {}))
        out["moments"] = {m: _Moments.merge_dicts([d["moments"][m] for d in dicts])
                          for m in ("plddt", "ptm")}
        rng = random.Random(seed)
        size = max(d["sample"]["size"] for d in dicts)
        pool = [(d["sample"]["rows"], d["sample"]["seen"]) for d in dicts if d["sample"]["rows"]]
        rows = []
        if pool:
            total_w = sum(w for _, w in pool)
            for _ in range(min(size, sum(len(p) for p, _ in pool))):
                r = rng.random() * total_w
                for p, w in pool:
                    r -= w
                    if r <= 0:
                        rows.append(p[rng.randrange(len(p))])
                        break
        out["hist"] = {k: v.tolist() for k, v in out["hist"].items()}
        out["status"] = dict(out["status"])
        out.update({k: dict(c) for k, c in counters.items()})
        out["sample"] = {"size": size, "seen": sum(w for _, w in pool), "rows": rows}
        return out


def quantiles_from_hist(counts, spec, qs=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)):
    """Quantiles read off a histogram -- resolution-limited to one bin width (0.5 pLDDT / 0.005 pTM),
    far finer than any figure resolves. Exact quantiles would need every value kept."""
    counts = np.asarray(counts, dtype=np.int64)
    lo, hi, nb = spec
    edges = np.linspace(lo, hi, nb + 1)
    tot = counts.sum()
    if not tot:
        return {str(q): None for q in qs}
    cum = np.cumsum(counts)
    out = {}
    for q in qs:
        k = min(int(np.searchsorted(cum, q * tot, side="left")), len(counts) - 1)
        out[str(q)] = float(edges[min(k, len(edges) - 1)])
    return out


# ---------------------------------------------------------------------------
# ledger (append-only, resumable)
# ---------------------------------------------------------------------------
COLS = ["orfid", "status", "seq_len", "truncated", "mean_plddt", "ptm", "iptm",
        "wall_s", "res_per_s"]


def read_done(path):
    """ORFids already in the ledger, so a restart resumes instead of redoing work."""
    if not os.path.exists(path):
        return set()
    done = set()
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("orfid"):
                done.add(row["orfid"])
    return done


def fetch_one(oid, src):
    rec = {c: "" for c in COLS}
    rec.update(orfid=oid, status="error")
    vals, status = src.scalars(oid)
    rec["status"] = status
    for k_out, k_in in (("seq_len", "seq_len"), ("truncated", "truncated"),
                        ("mean_plddt", "mean_plddt"), ("ptm", "ptm"), ("iptm", "iptm"),
                        ("wall_s", "wall_s"), ("res_per_s", "res_per_s")):
        v = vals.get(k_in)
        if v is not None:
            rec[k_out] = v
    return rec


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def build_worklist(a, src):
    """ORFs to read. The metrics listing is the DEFAULT and the authoritative source (an object
    exists iff the fold finished and uploaded); --ids only narrows it to a chosen subset."""
    if a.ids:
        ids = {norm_id(x) for x in open(a.ids).read().split() if x.strip()}
        return sorted(ids), f"--ids {os.path.basename(a.ids)}"
    return sorted(src.list_ids()), "metrics/ listing"


def run(a):
    shard = parse_shard(a.shard if a.shard is not None else os.environ.get("CONFMETRICS_SHARD"))
    tag = f".shard{shard[0]}of{shard[1]}" if shard else ""
    os.makedirs(a.out, exist_ok=True)
    ledger = os.path.join(a.out, f"conf{tag}.tsv")
    agg_path = os.path.join(a.out, f"aggregates{tag}.json")

    src = MetricSource(a.run, a.json_name, a.range_bytes, a.full_get)
    warm_imports(src.is_s3)                               # MUST precede the thread pool

    t_list = time.perf_counter()
    ids, how = build_worklist(a, src)
    print(f"run   = {a.label or src.label}   {a.run}")
    print(f"{len(ids)} folded ORFs (work list from {how}, {time.perf_counter() - t_list:.1f}s)")
    if shard:
        before = len(ids)
        ids = select_shard(ids, *shard)
        print(f"shard {shard[0]}/{shard[1]}: {len(ids)}/{before} ORFs")
    if a.limit:
        ids = ids[:a.limit]
        print(f"--limit {a.limit}: {len(ids)} ORFs")
    if not ids:
        raise SystemExit("nothing to read -- is the prefix right, and does it have a metrics/ dir?")

    done = read_done(ledger) if a.resume else set()
    if done:
        ids = [i for i in ids if i not in done]
        print(f"resume: {len(done)} already in {os.path.basename(ledger)}, {len(ids)} left")
    if not ids:
        print("nothing left to do -- ledger already complete for this shard")

    agg = Aggregates(a.sample, a.seed)
    lock = threading.Lock()
    t0 = time.perf_counter()
    fresh = not os.path.exists(ledger) or os.path.getsize(ledger) == 0
    with open(ledger, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        if fresh:
            w.writeheader()

        def work(oid):
            try:
                return fetch_one(oid, src)
            except Exception as e:                        # noqa: BLE001  one bad object must not kill the shard
                r = {c: "" for c in COLS}
                r.update(orfid=oid, status="error", iptm=f"{type(e).__name__}: {e}"[:120])
                return r

        n = 0
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            for rec in pool.map(work, ids):
                with lock:
                    w.writerow(rec)
                    agg.add(rec)
                    n += 1
                    if n % a.flush_every == 0:
                        fh.flush()
                        _write_agg(agg_path, agg, a, src, t0)
                        r = n / max(time.perf_counter() - t0, 1e-9)
                        print(f"  {n}/{len(ids)}  {r:.0f} ORFs/s  ok={agg.n}", flush=True)

    # An interrupted run leaves a valid ledger; on a resumed pass the aggregates are rebuilt from it
    # rather than trusting an in-memory object that only saw this pass's rows.
    if a.resume and done:
        agg = rebuild_aggregates(ledger, a.sample, a.seed)
    _write_agg(agg_path, agg, a, src, t0)
    summarise(agg, a, src, ledger, agg_path, time.perf_counter() - t0)
    if a.s3_out:
        mirror(a.s3_out, [ledger, agg_path])
    return 0


def rebuild_aggregates(ledger, sample, seed):
    """Recompute aggregates from a ledger TSV -- streaming, one row at a time (never loads the file)."""
    agg = Aggregates(sample, seed)
    with open(ledger, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rec = dict(row)
            try:
                rec["seq_len"] = int(rec.get("seq_len") or 0)
            except ValueError:
                rec["seq_len"] = 0
            agg.add(rec)
    return agg


def _write_agg(path, agg, a, src, t0):
    d = agg.to_dict()
    d["meta"] = {"run": a.run, "label": a.label or src.label, "shard": a.shard,
                 "range_bytes": None if a.full_get else a.range_bytes,
                 "full_get": bool(a.full_get), "widened": src.n_widened,
                 # 3 dp, not 1: a fast run (local dir, or a small shard) rounds to 0.0 at 1 dp and
                 # the notebook's throughput extrapolation then divides by zero.
                 "workers": a.workers, "elapsed_s": round(time.perf_counter() - t0, 3),
                 "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(d, f, separators=(",", ":"))
    os.replace(tmp, path)                                 # atomic -> a reader never sees half a file


def mirror(s3_out, paths):
    import boto3
    c = boto3.client("s3")
    bucket, _, prefix = s3_out[len("s3://"):].partition("/")
    for p in paths:
        key = "/".join(x for x in (prefix.strip("/"), os.path.basename(p)) if x)
        with open(p, "rb") as f:
            c.put_object(Bucket=bucket, Key=key, Body=f.read())
        print(f"  mirrored -> s3://{bucket}/{key}")


def summarise(agg, a, src, ledger, agg_path, elapsed):
    d = agg.to_dict()
    n = agg.n
    print("\n" + "=" * 72)
    print(f"Da  self-confidence  '{a.label or src.label}'    n_ok={n}   ({elapsed:.1f}s)")
    print("=" * 72)
    if n:
        for name, mom, hist, spec, thr, above, fmt in (
                ("mean pLDDT", d["moments"]["plddt"], d["hist"]["plddt"], PLDDT_BINS,
                 PLDDT_THRESHOLDS, agg.plddt_above, "{:.2f}"),
                ("pTM", d["moments"]["ptm"], d["hist"]["ptm"], PTM_BINS,
                 PTM_THRESHOLDS, agg.ptm_above, "{:.4f}")):
            q = quantiles_from_hist(hist, spec)
            print(f"\n  {name}")
            print(("    mean=" + fmt + "  sd=" + fmt + "  min=" + fmt + "  max=" + fmt).format(
                mom["mean"], math.sqrt(mom["var"]), mom["min"], mom["max"]))
            print("    quantiles: " + "  ".join(
                f"p{int(float(k) * 100)}=" + fmt.format(v) for k, v in q.items() if v is not None))
            print("    fraction at/above:  " + "  ".join(
                f"{t:g}={above[t]}/{n} ({100.0 * above[t] / n:.2f}%)" for t in thr))
        if agg.n_truncated:
            t70 = agg.trunc_plddt_above[70.0]
            print(f"\n  {agg.n_truncated} ORFs were TRUNCATED to the ESMC context before folding "
                  f"({100.0 * agg.n_truncated / n:.2f}%);")
            print(f"    their confidence describes the folded prefix, not the ORF "
                  f"({t70}/{agg.n_truncated} at/above pLDDT 70).")
    bad = {k: v for k, v in d["status"].items() if k != "ok"}
    if bad:
        print("\n  non-ok: " + ", ".join(f"{k}={v}" for k, v in sorted(bad.items())))
    if src.n_widened:
        print(f"  {src.n_widened} object(s) needed a full GET (head window too small) "
              f"-- raise --range-bytes if this is not ~0")
    print("\n  NOTE the denominator is structures that EXIST. Folds that errored never wrote a")
    print("  metrics JSON, so the attempt-level failure rate is not visible here (by design).")
    print(f"\n  per-ORF rows -> {ledger}")
    print(f"  aggregates   -> {agg_path}   (plot from this; O(1) in N)")
    print("=" * 72)


def merge_mode(a):
    """Combine per-shard aggregates*.json into one aggregates.merged.json for plotting."""
    import glob
    paths = sorted(glob.glob(os.path.join(a.out, "aggregates.shard*.json")))
    if not paths:
        raise SystemExit(f"no aggregates.shard*.json in {a.out}")
    merged = Aggregates.merge([json.load(open(p)) for p in paths], a.seed)
    merged["meta"] = {"merged_from": [os.path.basename(p) for p in paths],
                      "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    out = os.path.join(a.out, "aggregates.merged.json")
    with open(out, "w") as f:
        json.dump(merged, f, separators=(",", ":"))
    print(f"merged {len(paths)} shard aggregates  n_ok={merged['n_ok']}  -> {out}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="prediction run: s3://bucket/prefix/ or a local dir "
                                  "(reads <root>/metrics/*.json)")
    ap.add_argument("--out", default="conf", help="output dir for ledger + aggregates (default conf/)")
    ap.add_argument("--label", default="", help="short name for this run in plots/summaries")
    ap.add_argument("--ids", help="file of ORFids to read (one per token; 'orf' prefix optional). "
                                  "Default is the metrics/ listing, which is authoritative.")
    ap.add_argument("--json-name", default="orf{id}.json", dest="json_name",
                    help="metrics filename template (default orf{id}.json)")
    ap.add_argument("--range-bytes", type=int, default=4096,
                    help="bytes of each object to fetch; the scalars sit at the front, before the "
                         "O(L^2) pae array (default 4096). Objects whose head lacks a required "
                         "scalar are automatically re-fetched whole.")
    ap.add_argument("--full-get", action="store_true",
                    help="fetch whole objects instead of the head (~500x the bytes; only needed if "
                         "the metrics layout changes)")
    ap.add_argument("--workers", type=int, default=32,
                    help="concurrent GETs (default 32). This run is request-latency-bound, so more "
                         "workers scale it nearly linearly -- push to 64-128 on an in-region box.")
    ap.add_argument("--shard", default=None, metavar="K/N",
                    help="process only band K of N (1-based) -- one command, a fleet of workers. "
                         "Also read from $CONFMETRICS_SHARD.")
    ap.add_argument("--limit", type=int, help="only read the first N ORFs (smoke test)")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                    help="ignore an existing ledger and re-read everything")
    ap.add_argument("--sample", type=int, default=20000,
                    help="reservoir subsample kept in the aggregates for scatter panels (default 20k)")
    ap.add_argument("--seed", type=int, default=0, help="reservoir seed (reproducible subsample)")
    ap.add_argument("--flush-every", type=int, default=1000,
                    help="rows between ledger flush + aggregates rewrite (default 1000)")
    ap.add_argument("--s3-out", metavar="s3://bucket/prefix/",
                    help="also mirror the ledger + aggregates here when the run finishes")
    ap.add_argument("--merge", action="store_true",
                    help="merge aggregates.shard*.json in --out into aggregates.merged.json and exit")
    a = ap.parse_args()
    if a.merge:
        return merge_mode(a)
    if not a.run:
        ap.error("--run is required (or use --merge)")
    return run(a)


if __name__ == "__main__":
    sys.exit(main())
