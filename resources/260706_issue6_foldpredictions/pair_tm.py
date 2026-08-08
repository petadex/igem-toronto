"""
pair_tm.py -- figure Ba: TM-score between TWO predictions of the SAME protein.

Given two prediction runs that folded the same ORFs (run A = baseline, run B = the variant, e.g.
MSA-conditioned), score every shared ORF pair and answer "did the change move the fold at all?".
This is a run-vs-run comparison -- NO ground-truth PDB is involved, so it covers every folded ORF,
not just the handful with a solved structure.

    python pair_tm.py --run-a s3://bucket/base/ --run-b s3://bucket/msa/ --out pairtm_base_vs_msa

Both --run-a/--run-b are either an `s3://bucket/prefix/` or a local directory, laid out the way
esmfold2_local_predictor.py writes them:  <root>/structures/orf<ORFid>.cif  (+ tracker.csv).

BUILT FOR THE FULL CENTROID SET (millions of pairs), so nothing here scales with N:
  * structures are STREAMED (one GET each, parsed from memory -- never written to disk);
  * the work list is an explicit id list (or the two runs' tracker.csv), not an S3 listing walk;
  * `--shard K/N` fans one command across a fleet; each shard owns a disjoint id band and its own
    ledger file, so workers never clobber each other (same convention as the predictor's --shard);
  * the ledger is append-only and re-read on start, so an interrupted run RESUMES instead of redoing;
  * per-ORF rows are a convenience, not the plotting input -- the real output is `aggregates*.json`
    (fixed-bin histograms + Welford moments + exact threshold counts + a reservoir subsample). Those
    are O(1) in N and MERGEABLE across shards (`--merge`), so the figures cost the same at 121 pairs
    or 2 million.

READ THIS BEFORE INTERPRETING:
  * TM here measures how much the prediction MOVED between the two runs, not whether it got better.
    TM ~= 1 means the change was inert; a low TM means "different fold", not "worse fold". Accuracy
    needs experimental references (that is group C, a different script).
  * Both normalisations are reported (tm_norm_a / tm_norm_b). For same-sequence pairs the two chains
    are the same length so they agree; they diverge only if a run truncated differently.

CHOICE OF METRIC (--metric auto, the default):
  Two folds of the SAME sequence have a KNOWN 1:1 residue correspondence, so a Kabsch superposition
  over that correspondence is the exact answer -- not an approximation. TM-align's freedom to find its
  own correspondence is only *needed* when the sequences differ (a run truncated, or the wrong prefix
  was passed). So `auto` uses kabsch when the two sequences match and falls back to tm-align when they
  do not, recording which it used per row in `metric_used`.
  This matters because the two differ enormously in cost -- measured on this pipeline's structures:

      length     tm-align      kabsch
       150 aa      12 ms       0.12 ms
       300 aa      40 ms       0.10 ms
       600 aa     330 ms       0.11 ms
       900 aa    1351 ms       0.11 ms

  TM-align is ~cubic in length while kabsch is flat, so the ORFs near the 2048 aa ESMC context would
  dominate a full-set run on their own. On real predictions the two agree to ~5e-3 (see --metric both
  and the notebook's proxy cell). NOTE kabsch is NOT a strict lower bound on tm-align, despite having
  a fixed correspondence: tm-align's superposition is an iterative heuristic, not a global optimum, so
  it can land marginally BELOW the kabsch value on near-identical structures. Treat the difference as
  two-sided agreement of order 1e-3, not as a bound.

Deps:  pip install tmtools numpy boto3      (biopython only for --parser biopython)
"""
import argparse
import csv
import io
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
# mmCIF -> (CA coords, one-letter sequence)
# ---------------------------------------------------------------------------
# Our predictions come from our own writer, so the default parser is a plain tokenizer over the
# _atom_site loop (~1 ms/structure vs ~20 ms for Biopython). It still reads the column order from
# the loop header rather than assuming it, so it also handles real mmCIF. --parser biopython swaps in
# the strict reader; the self-test asserts the two agree.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# Anything else -> "X" (never dropped). ESMFold2 emits UNK wherever the ORF had an ambiguous X;
# dropping those residues would silently truncate the chain, and TM-align is coordinate-based, so
# an occasional X is harmless for the superposition.


def parse_cif_fast(text):
    """(coords[N,3] float64, seq str) for CA atoms in label_seq_id order. ('', empty) if none."""
    lines = text.splitlines()
    i, n = 0, len(lines)
    coords, resn, seqid = [], [], []
    while i < n:
        if lines[i].strip() != "loop_":
            i += 1
            continue
        cols, j = [], i + 1
        while j < n and lines[j].lstrip().startswith("_"):
            cols.append(lines[j].split()[0].strip())    # header token only (ignore trailing spaces)
            j += 1
        if not any(c.startswith("_atom_site.") for c in cols):
            i = j
            continue
        idx = {c: k for k, c in enumerate(cols)}
        need = ["_atom_site.label_atom_id", "_atom_site.label_comp_id",
                "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
        if not all(c in idx for c in need):
            i = j
            continue
        i_at, i_cp = idx["_atom_site.label_atom_id"], idx["_atom_site.label_comp_id"]
        i_x, i_y, i_z = (idx["_atom_site.Cartn_x"], idx["_atom_site.Cartn_y"], idx["_atom_site.Cartn_z"])
        i_sq = idx.get("_atom_site.label_seq_id")
        ncol = len(cols)
        k = j
        while k < n:
            s = lines[k].strip()
            if not s or s[0] in "#_" or s.startswith(("loop_", "data_")):
                break
            tok = s.split()
            if len(tok) >= ncol and tok[i_at].strip('"') == "CA":
                try:
                    xyz = (float(tok[i_x]), float(tok[i_y]), float(tok[i_z]))
                except ValueError:
                    k += 1
                    continue
                coords.append(xyz)
                resn.append(tok[i_cp].strip('"').upper())
                seqid.append(_asint(tok[i_sq], len(coords)) if i_sq is not None else len(coords))
            k += 1
        i = k
    if not coords:
        return np.zeros((0, 3)), ""
    order = np.argsort(np.asarray(seqid), kind="stable")
    coords = np.asarray(coords, dtype=float)[order]
    seq = "".join(THREE_TO_ONE.get(resn[o], "X") for o in order)
    return coords, seq


def _asint(x, fallback):
    try:
        return int(x)
    except (ValueError, TypeError):
        return fallback


def parse_cif_biopython(text):
    """Strict-reader cross-check. Same contract as parse_cif_fast; parses from memory (no temp file).

    NOT the safer option, despite being the "real" parser: it requires _atom_site.id and
    _atom_site.label_alt_id and KeyErrors without them, so it rejects ESMAtlas-style minimal CIFs
    that parse_cif_fast reads fine (verified on atlastest/*.cif). It exists to validate the fast
    parser on OUR predictions -- where the two agree exactly, at ~240x the cost.
    """
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.SeqUtils import seq1
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")                 # our minimal CIFs lack header blocks
        s = MMCIFParser(QUIET=True).get_structure("p", io.StringIO(text))
    model = next(s.get_models())
    best = (np.zeros((0, 3)), "")
    for chain in model.get_chains():                    # predictions are single-chain; keep longest
        co, sq = [], []
        for res in chain.get_residues():
            if "CA" in res:
                co.append(res["CA"].coord)
                sq.append(seq1(res.resname, undef_code="X") or "X")
        if len(sq) > len(best[1]):
            best = (np.asarray(co, dtype=float), "".join(sq))
    return best


PARSERS = {"fast": parse_cif_fast, "biopython": parse_cif_biopython}


def warm_imports(metric, parser, needs_s3):
    """Import every deferred dependency ONCE, in the main thread, BEFORE the worker pool starts.

    This is load-bearing, not tidiness. tmtools / Bio / boto3 are imported lazily (so the script runs
    with only the deps a given mode needs), but if the FIRST import happens concurrently in several
    worker threads they deadlock on the import lock -- the whole run hangs with every thread parked
    inside `_find_and_load`. Warming them here makes the in-function imports a sys.modules hit.

    tmtools needs a dummy CALL, not just an import: `tm_align` pulls in its compiled extension on
    first invocation, so importing it in the main thread but first CALLING it from N threads
    deadlocks exactly the same way (this is what actually hung -- the import frames sit *inside*
    tm_align, above our call site).
    """
    if metric in ("auto", "tmalign", "both"):            # 'auto' may still need the fallback
        from tmtools import tm_align
        global _TM_ALIGN
        _TM_ALIGN = tm_align
        c = np.stack([np.zeros(10), np.zeros(10), np.arange(10) * 3.8], 1)
        tm_align(c, c + 1.0, "A" * 10, "A" * 10)         # force tmtools' own lazy import
    if parser == "biopython":
        from Bio.PDB.MMCIFParser import MMCIFParser      # noqa: F401
        from Bio.SeqUtils import seq1                    # noqa: F401
        parse_cif_biopython(_fake_cif(np.zeros((3, 3)), "AAA"))   # force Bio's lazy machinery
    if needs_s3:
        import boto3
        boto3.client("s3")                               # builds botocore's loaders/session once


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
_TM_ALIGN = None      # resolved once by warm_imports(); see the deadlock note there


def tm_align_pair(ca, sa, cb, sb):
    """TM-align: (tm_norm_a, tm_norm_b, rmsd). Sequence-independent, finds its own correspondence."""
    global _TM_ALIGN
    if _TM_ALIGN is None:
        from tmtools import tm_align
        _TM_ALIGN = tm_align
    r = _TM_ALIGN(ca, cb, sa, sb)
    return float(r.tm_norm_chain1), float(r.tm_norm_chain2), float(r.rmsd)


def kabsch_rmsd_tm(P, Q):
    """RMSD + TM-score over the KNOWN 1:1 correspondence (equal-length inputs). Rotation/translation
    invariant, O(L), and EXACT when the correspondence really is known -- which it is for two folds of
    the same sequence. See the module docstring on why this is not a strict bound on TM-align."""
    L = len(P)
    if L == 0:
        return float("nan"), float("nan")
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    di = np.sqrt((((Pc @ R.T) - Qc) ** 2).sum(1))
    d0 = max(1.24 * (L - 15) ** (1.0 / 3.0) - 1.8, 0.5) if L > 15 else 0.5
    return float(np.sqrt((di ** 2).mean())), float(np.mean(1.0 / (1.0 + (di / d0) ** 2)))


# ---------------------------------------------------------------------------
# structure source: S3 stream (default) or local dir
# ---------------------------------------------------------------------------
class StructSource:
    """Reads <root>/structures/<cif_name> as text. `root` is s3://bucket/prefix/ or a local dir.

    One boto3 client PER THREAD: clients are only nominally thread-safe and sharing one across a
    large pool serialises on its connection pool.
    """

    def __init__(self, root, cif_name="orf{id}.cif", retries=4):
        self.root, self.cif_name, self.retries = root, cif_name, retries
        self.is_s3 = str(root).startswith("s3://")
        self._tl = threading.local()
        if self.is_s3:
            self.bucket, _, self.prefix = root[len("s3://"):].partition("/")
            self.prefix = self.prefix.strip("/")
        self.label = os.path.basename(str(root).rstrip("/\\")) or str(root)

    @property
    def s3(self):
        if getattr(self._tl, "c", None) is None:
            import boto3                                 # pre-warmed by warm_imports()
            self._tl.c = boto3.client("s3")
        return self._tl.c

    def _key(self, oid):
        return "/".join(p for p in (self.prefix, "structures", self.cif_name.format(id=oid)) if p)

    def get_text(self, oid):
        """CIF text, or None if this run has no structure for `oid` (not folded / not uploaded yet)."""
        name = self.cif_name.format(id=oid)
        if not self.is_s3:
            for p in (os.path.join(self.root, "structures", name), os.path.join(self.root, name)):
                if os.path.exists(p):
                    with open(p) as f:
                        return f.read()
            return None
        for attempt in range(self.retries):
            try:
                return self.s3.get_object(Bucket=self.bucket, Key=self._key(oid))["Body"].read().decode()
            except Exception as e:                       # noqa: BLE001
                code = str(getattr(e, "response", {}).get("Error", {}).get("Code", ""))
                if code in ("NoSuchKey", "404", "NoSuchBucket") or type(e).__name__ == "NoSuchKey":
                    return None
                if attempt == self.retries - 1:
                    raise
                time.sleep(0.4 * 2 ** attempt)           # throttling/transient -> backoff
        return None

    def list_ids(self):
        """Every ORFid with a structure under this root (paginated; 1000 keys/request)."""
        pat = re.compile("^" + re.escape(self.cif_name).replace(r"\{id\}", "(.+)") + "$")
        out = set()
        if not self.is_s3:
            d = os.path.join(self.root, "structures")
            for n in (os.listdir(d) if os.path.isdir(d) else []):
                m = pat.match(n)
                if m:
                    out.add(m.group(1))
            return out
        tok = None
        while True:
            kw = {"Bucket": self.bucket, "Prefix": "/".join(p for p in (self.prefix, "structures/") if p)}
            if tok:
                kw["ContinuationToken"] = tok
            r = self.s3.list_objects_v2(**kw)
            for o in r.get("Contents", []):
                m = pat.match(o["Key"].rsplit("/", 1)[-1])
                if m:
                    out.add(m.group(1))
            if not r.get("IsTruncated"):
                return out
            tok = r["NextContinuationToken"]

    def tracker_ids(self):
        """ORFids from this run's tracker.csv (or merged tracker*.csv shards). Empty set if absent."""
        texts = []
        if self.is_s3:
            names = []
            tok = None
            while True:
                kw = {"Bucket": self.bucket, "Prefix": "/".join(p for p in (self.prefix, "tracker") if p)}
                if tok:
                    kw["ContinuationToken"] = tok
                r = self.s3.list_objects_v2(**kw)
                names += [o["Key"] for o in r.get("Contents", []) if o["Key"].endswith(".csv")]
                if not r.get("IsTruncated"):
                    break
                tok = r["NextContinuationToken"]
            for k in names:
                texts.append(self.s3.get_object(Bucket=self.bucket, Key=k)["Body"].read().decode())
        else:
            import glob
            for p in sorted(glob.glob(os.path.join(self.root, "tracker*.csv"))):
                with open(p) as f:
                    texts.append(f.read())
        out = set()
        for t in texts:
            for row in csv.DictReader(io.StringIO(t)):
                if row.get("id"):
                    out.add(norm_id(row["id"]))
        return out


def norm_id(x):
    """Canonical bare ORFid. tracker.csv ids look like 'orf3772973' while the map files use bare
    '3772973'; both must land on the same key so ledgers and joins line up."""
    x = str(x).strip()
    return x[3:] if x.startswith("orf") else x


# ---------------------------------------------------------------------------
# sharding (same 'K/N' contract as esmfold2_local_predictor.py --shard)
# ---------------------------------------------------------------------------
def parse_shard(spec):
    if not spec:
        return None
    try:
        ks, ns = str(spec).split("/")
        k, n = int(ks), int(ns)
    except Exception:                                    # noqa: BLE001
        raise SystemExit(f"--shard must be 'K/N' (1-based), got {spec!r}")
    if not (1 <= k <= n):
        raise SystemExit(f"--shard 'K/N' needs 1<=K<=N, got {k}/{n}")
    return k, n


def select_shard(ids, k, n):
    """Contiguous band K of N over ids sorted deterministically -> shards never overlap or gap,
    whatever order the ids arrived in. Unlike the predictor we do NOT length-band (no compile cache
    to keep warm here); numeric-aware id order is enough and keeps bands reproducible."""
    ordered = sorted(ids, key=lambda s: (len(s), s))
    b = [round(i * len(ordered) / n) for i in range(n + 1)]
    return ordered[b[k - 1]:b[k]]


# ---------------------------------------------------------------------------
# streaming aggregates -- the actual plotting input. O(1) in N, mergeable across shards.
# ---------------------------------------------------------------------------
# Histogram edges are FIXED (never data-derived) so shard JSONs stay bin-for-bin addable.
TM_BINS = (0.0, 1.0, 200)          # 0.005 wide
TM_ZOOM_BINS = (0.9, 1.0, 200)     # 0.0005 wide -- run-vs-run TM piles up near 1.0, and 6a showed
                                   # the interesting spread lives inside the last percent
RMSD_BINS = (0.0, 10.0, 200)       # 0.05 A wide, + overflow bucket
LEN_BINS = (0.0, 2100.0, 42)       # 50 aa wide (ESMC context is 2048)
THRESHOLDS = (0.5, 0.7, 0.8, 0.9, 0.95, 0.99)


def _edges(spec):
    lo, hi, nb = spec
    return np.linspace(lo, hi, nb + 1)


class Aggregates:
    """Fixed-bin histograms + Welford moments + exact threshold counts + a reservoir subsample."""

    def __init__(self, sample_size=20000, seed=0):
        self.n = 0
        self.tm = np.zeros(TM_BINS[2] + 1, dtype=np.int64)        # last cell = >= hi (i.e. tm == 1)
        self.tm_zoom = np.zeros(TM_ZOOM_BINS[2] + 2, dtype=np.int64)  # [0] = < 0.9, [-1] = >= 1.0
        self.rmsd = np.zeros(RMSD_BINS[2] + 1, dtype=np.int64)    # last cell = overflow
        self.tm_by_len = np.zeros((LEN_BINS[2] + 1, 100), dtype=np.int64)
        self.below = {t: 0 for t in THRESHOLDS}                   # exact, not read off the bins
        self.mean = 0.0
        self._m2 = 0.0
        self.min = float("inf")
        self.max = float("-inf")
        self.status = Counter()
        self.sample_size = sample_size
        self.sample = []
        self._seen = 0                                            # reservoir stream position
        self._rng = random.Random(seed)

    # -- one finished pair -------------------------------------------------
    def add(self, rec):
        self.status[rec["status"]] += 1
        if rec["status"] != "ok":
            return
        tm, rmsd, L = rec["tm_min"], rec["rmsd"], rec["n_res_a"]
        self.n += 1
        self.tm[_bin(tm, TM_BINS)] += 1
        self.tm_zoom[0 if tm < 0.9 else (_bin(tm, TM_ZOOM_BINS) + 1)] += 1
        self.rmsd[_bin(rmsd, RMSD_BINS)] += 1
        self.tm_by_len[_bin(float(L), LEN_BINS), min(int(tm * 100), 99)] += 1
        for t in THRESHOLDS:
            if tm < t:
                self.below[t] += 1
        d = tm - self.mean                                        # Welford
        self.mean += d / self.n
        self._m2 += d * (tm - self.mean)
        self.min, self.max = min(self.min, tm), max(self.max, tm)
        self._reservoir(rec)

    def _reservoir(self, rec):
        """Algorithm R, seeded -> a uniform sample of the stream for scatter/outlier inspection at
        any N, reproducible across reruns of the same shard."""
        self._seen += 1
        row = {k: rec[k] for k in ("orfid", "tm_min", "rmsd", "n_res_a")}
        if len(self.sample) < self.sample_size:
            self.sample.append(row)
        else:
            j = self._rng.randrange(self._seen)
            if j < self.sample_size:
                self.sample[j] = row

    # -- serialise / merge -------------------------------------------------
    def to_dict(self):
        return {
            "n_ok": self.n, "status": dict(self.status),
            "bins": {"tm": TM_BINS, "tm_zoom": TM_ZOOM_BINS, "rmsd": RMSD_BINS, "len": LEN_BINS},
            "hist": {"tm": self.tm.tolist(), "tm_zoom": self.tm_zoom.tolist(),
                     "rmsd": self.rmsd.tolist(), "tm_by_len": self.tm_by_len.tolist()},
            "tm_below": {str(k): v for k, v in self.below.items()},
            "moments": {"mean": self.mean, "m2": self._m2,
                        "var": (self._m2 / (self.n - 1)) if self.n > 1 else 0.0,
                        "min": self.min if self.n else None, "max": self.max if self.n else None},
            "sample": {"size": self.sample_size, "seen": self._seen, "rows": self.sample},
        }

    @staticmethod
    def merge(dicts, seed=0):
        """Combine shard aggregates into one. Histograms add; Welford combines by the parallel
        formula; reservoirs are resampled from the union weighted by each shard's stream length
        (each shard's sample is uniform over its own stream, so this is uniform over the union)."""
        dicts = [d for d in dicts if d]
        if not dicts:
            raise SystemExit("nothing to merge")
        out = {"n_ok": 0, "status": Counter(), "bins": dicts[0]["bins"],
               "hist": {k: np.zeros_like(np.asarray(v, dtype=np.int64))
                        for k, v in dicts[0]["hist"].items()},
               "tm_below": Counter()}
        mean, m2, n = 0.0, 0.0, 0
        mn, mx = float("inf"), float("-inf")
        for d in dicts:
            if d["bins"] != out["bins"]:
                raise SystemExit("bin edges differ between shards -- cannot merge")
            out["n_ok"] += d["n_ok"]
            out["status"].update(d["status"])
            for k in out["hist"]:
                out["hist"][k] += np.asarray(d["hist"][k], dtype=np.int64)
            out["tm_below"].update({k: v for k, v in d["tm_below"].items()})
            dn, dmean, dm2 = d["n_ok"], d["moments"]["mean"], d["moments"]["m2"]
            if dn:
                delta = dmean - mean
                tot = n + dn
                mean += delta * dn / tot
                m2 += dm2 + delta ** 2 * n * dn / tot
                n = tot
                mn = min(mn, d["moments"]["min"])
                mx = max(mx, d["moments"]["max"])
        rng = random.Random(seed)
        size = max(d["sample"]["size"] for d in dicts)
        pool, weights = [], []
        for d in dicts:
            rows = d["sample"]["rows"]
            if rows:
                pool.append(rows)
                weights.append(d["sample"]["seen"])
        rows = []
        if pool:
            total_w = sum(weights)
            for _ in range(min(size, sum(len(p) for p in pool))):
                r = rng.random() * total_w
                for p, w in zip(pool, weights):
                    r -= w
                    if r <= 0:
                        rows.append(p[rng.randrange(len(p))])
                        break
        out["hist"] = {k: v.tolist() for k, v in out["hist"].items()}
        out["status"] = dict(out["status"])
        out["tm_below"] = dict(out["tm_below"])
        out["moments"] = {"mean": mean, "m2": m2, "var": (m2 / (n - 1)) if n > 1 else 0.0,
                          "min": mn if n else None, "max": mx if n else None}
        out["sample"] = {"size": size, "seen": sum(weights), "rows": rows}
        return out


def _bin(v, spec):
    """Index into `spec`'s bins, with everything >= hi landing in one extra overflow cell."""
    lo, hi, nb = spec
    if not (v == v):                                     # NaN
        return nb
    if v >= hi:
        return nb
    return max(0, min(nb - 1, int((v - lo) / (hi - lo) * nb)))


def quantiles_from_hist(counts, spec, qs=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), zoom_lo=None):
    """Quantiles read off a histogram -- resolution-limited to one bin width, which is 0.005 (or
    0.0005 on the zoom), far finer than any figure resolves. Exact quantiles would need every value."""
    counts = np.asarray(counts, dtype=np.int64)
    lo, hi, nb = spec
    edges = np.linspace(lo, hi, nb + 1)
    if zoom_lo is not None:                              # zoom hist carries an underflow cell at [0]
        counts, edges = counts[1:], np.append(edges, hi)
    tot = counts.sum()
    if not tot:
        return {str(q): None for q in qs}
    cum = np.cumsum(counts)
    out = {}
    for q in qs:
        k = int(np.searchsorted(cum, q * tot, side="left"))
        k = min(k, len(counts) - 1)
        out[str(q)] = float(edges[min(k, len(edges) - 1)])
    return out


# ---------------------------------------------------------------------------
# ledger (append-only, resumable)
# ---------------------------------------------------------------------------
COLS = ["orfid", "status", "metric_used", "n_res_a", "n_res_b", "seq_match",
        "tm_norm_a", "tm_norm_b", "tm_min", "rmsd", "tm_kabsch", "rmsd_kabsch"]


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


# ---------------------------------------------------------------------------
# scoring one pair
# ---------------------------------------------------------------------------
def score_pair(oid, src_a, src_b, parse, metric):
    rec = {c: "" for c in COLS}
    rec.update(orfid=oid, status="error", tm_min=float("nan"), rmsd=float("nan"),
               n_res_a=0, n_res_b=0)

    ta = src_a.get_text(oid)
    if ta is None:
        rec["status"] = "missing_a"
        return rec
    tb = src_b.get_text(oid)
    if tb is None:
        rec["status"] = "missing_b"
        return rec

    ca, sa = parse(ta)
    cb, sb = parse(tb)
    if len(ca) == 0:
        rec["status"] = "parse_fail_a"
        return rec
    if len(cb) == 0:
        rec["status"] = "parse_fail_b"
        return rec
    rec["n_res_a"], rec["n_res_b"] = len(ca), len(cb)
    # Same ORF folded twice -> the sequences MUST match. A mismatch means the two runs were not fed
    # the same input (wrong prefix, stale fold, different truncation); flag it rather than quietly
    # reporting a TM between two different proteins.
    rec["seq_match"] = "exact" if sa == sb else ("len_diff" if len(sa) != len(sb) else "seq_diff")

    # 'auto': the correspondence is KNOWN when the sequences match, so kabsch is exact and ~1000x
    # cheaper; only fall back to TM-align when it isn't (truncation / mismatched inputs).
    if metric == "auto":
        metric = "kabsch" if rec["seq_match"] == "exact" else "tmalign"
    rec["metric_used"] = metric

    if metric in ("tmalign", "both"):
        try:
            tm_a, tm_b, rmsd = tm_align_pair(ca, sa, cb, sb)
        except Exception as e:                           # noqa: BLE001
            rec["status"] = "align_fail"
            rec["tm_kabsch"] = f"{type(e).__name__}"
            return rec
        rec.update(tm_norm_a=round(tm_a, 5), tm_norm_b=round(tm_b, 5),
                   tm_min=min(tm_a, tm_b), rmsd=round(rmsd, 4))
    if metric in ("kabsch", "both"):
        if len(ca) == len(cb):
            rk, tk = kabsch_rmsd_tm(ca, cb)
            rec["tm_kabsch"], rec["rmsd_kabsch"] = round(tk, 5), round(rk, 4)
            if metric == "kabsch":
                rec.update(tm_norm_a=round(tk, 5), tm_norm_b=round(tk, 5), tm_min=tk, rmsd=round(rk, 4))
        elif metric == "kabsch":
            rec["status"] = "len_mismatch_kabsch"        # no 1:1 correspondence -> use tmalign
            return rec

    rec["status"] = "ok"
    return rec


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def build_worklist(a, args, src_a, src_b):
    """Pairs to score = ids present in BOTH runs. Prefer an explicit list (cheapest, deterministic),
    then the runs' trackers, then an S3 listing walk."""
    if a.ids:
        ids = {norm_id(x) for x in open(a.ids).read().split() if x.strip()}
        src = f"--ids {a.ids}"
        if not a.no_intersect:
            pass                                         # trust the list; missing -> missing_a/b rows
        return sorted(ids), src, None, None
    if not a.ids_from_listing:
        ia, ib = src_a.tracker_ids(), src_b.tracker_ids()
        if ia and ib:
            return sorted(ia & ib), "tracker.csv of both runs", len(ia - ib), len(ib - ia)
    ia, ib = src_a.list_ids(), src_b.list_ids()
    ia, ib = {norm_id(x) for x in ia}, {norm_id(x) for x in ib}
    return sorted(ia & ib), "S3/dir listing of both runs", len(ia - ib), len(ib - ia)


def run(a):
    shard = parse_shard(a.shard if a.shard is not None else os.environ.get("PAIRTM_SHARD"))
    tag = f".shard{shard[0]}of{shard[1]}" if shard else ""
    os.makedirs(a.out, exist_ok=True)
    ledger = os.path.join(a.out, f"pairs{tag}.tsv")
    agg_path = os.path.join(a.out, f"aggregates{tag}.json")

    src_a = StructSource(a.run_a, a.cif_name)
    src_b = StructSource(a.run_b, a.cif_name)
    parse = PARSERS[a.parser]
    warm_imports(a.metric, a.parser, src_a.is_s3 or src_b.is_s3)   # MUST precede the thread pool

    ids, how, only_a, only_b = build_worklist(a, a, src_a, src_b)
    print(f"run A = {a.label_a or src_a.label}   {a.run_a}")
    print(f"run B = {a.label_b or src_b.label}   {a.run_b}")
    print(f"{len(ids)} shared ORFs (work list from {how})"
          + (f" | only in A: {only_a} | only in B: {only_b}" if only_a or only_b else ""))
    if shard:
        before = len(ids)
        ids = select_shard(ids, *shard)
        print(f"shard {shard[0]}/{shard[1]}: {len(ids)}/{before} pairs")
    if a.limit:
        ids = ids[:a.limit]
        print(f"--limit {a.limit}: {len(ids)} pairs")
    if not ids:
        raise SystemExit("no pairs to score")

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
                return score_pair(oid, src_a, src_b, parse, a.metric)
            except Exception as e:                       # noqa: BLE001  one bad pair must not kill the shard
                r = {c: "" for c in COLS}
                r.update(orfid=oid, status="error", tm_min=float("nan"), rmsd=float("nan"),
                         n_res_a=0, n_res_b=0, tm_kabsch=f"{type(e).__name__}: {e}"[:120])
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
                        _write_agg(agg_path, agg, a, src_a, src_b, t0)
                        r = n / max(time.perf_counter() - t0, 1e-9)
                        print(f"  {n}/{len(ids)}  {r:.1f} pairs/s  ok={agg.n}", flush=True)

    # An interrupted run leaves a valid ledger; the aggregates are rebuilt from it rather than
    # trusting a partial in-memory object (and to pick up rows from earlier resumed passes).
    if a.resume and done:
        agg = rebuild_aggregates(ledger, a.sample, a.seed)
    _write_agg(agg_path, agg, a, src_a, src_b, t0)
    summarise(agg, a, ledger, agg_path, time.perf_counter() - t0)
    if a.s3_out:
        mirror(a.s3_out, [ledger, agg_path])
    return 0


def rebuild_aggregates(ledger, sample, seed):
    """Recompute aggregates from a ledger TSV -- streaming, one row at a time (never loads the file)."""
    agg = Aggregates(sample, seed)
    with open(ledger, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rec = dict(row)
            for k in ("tm_min", "rmsd"):
                try:
                    rec[k] = float(rec[k])
                except (ValueError, TypeError):
                    rec[k] = float("nan")
            rec["n_res_a"] = _asint(rec.get("n_res_a"), 0)
            agg.add(rec)
    return agg


def _write_agg(path, agg, a, src_a, src_b, t0):
    d = agg.to_dict()
    d["meta"] = {"run_a": a.run_a, "run_b": a.run_b,
                 "label_a": a.label_a or src_a.label, "label_b": a.label_b or src_b.label,
                 "metric": a.metric, "parser": a.parser, "shard": a.shard,
                 "elapsed_s": round(time.perf_counter() - t0, 1),
                 "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(d, f, separators=(",", ":"))
    os.replace(tmp, path)                                # atomic -> a reader never sees half a file


def mirror(s3_out, paths):
    import boto3
    c = boto3.client("s3")
    bucket, _, prefix = s3_out[len("s3://"):].partition("/")
    for p in paths:
        key = "/".join(x for x in (prefix.strip("/"), os.path.basename(p)) if x)
        with open(p, "rb") as f:
            c.put_object(Bucket=bucket, Key=key, Body=f.read())
        print(f"  mirrored -> s3://{bucket}/{key}")


def summarise(agg, a, ledger, agg_path, elapsed):
    d = agg.to_dict()
    q = quantiles_from_hist(d["hist"]["tm"], TM_BINS)
    la, lb = a.label_a or "A", a.label_b or "B"
    print("\n" + "=" * 68)
    print(f"Ba  TM-score  '{la}'  vs  '{lb}'    n_ok={agg.n}   ({elapsed:.1f}s)")
    print("=" * 68)
    if agg.n:
        print(f"  mean={d['moments']['mean']:.4f}  sd={math.sqrt(d['moments']['var']):.4f}  "
              f"min={d['moments']['min']:.4f}  max={d['moments']['max']:.4f}")
        print("  quantiles (bin-resolution 0.005): "
              + "  ".join(f"p{int(float(k)*100)}={v:.3f}" for k, v in q.items() if v is not None))
        print("  fraction below:  " + "  ".join(
            f"{t}={agg.below[t]}/{agg.n} ({100.0*agg.below[t]/agg.n:.2f}%)" for t in THRESHOLDS))
    bad = {k: v for k, v in d["status"].items() if k != "ok"}
    if bad:
        print("  non-ok: " + ", ".join(f"{k}={v}" for k, v in sorted(bad.items())))
    if agg.n:
        print(f"\n  TM ~= 1 means the change was INERT (fold did not move), not that it was good.")
    print(f"\n  per-pair rows -> {ledger}")
    print(f"  aggregates    -> {agg_path}   (plot from this; O(1) in N)")
    print("=" * 68)


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


# ---------------------------------------------------------------------------
# self-test -- fabricates two runs on disk; no S3, no creds, no GPU
# ---------------------------------------------------------------------------
def _fake_cif(coords, seq):
    """A CA-only mmCIF using esmfold2_local_predictor.py's EXACT _atom_site column set and order
    (Cartn_* late, id last) -- so the test exercises header-driven column lookup on the real layout,
    not a convenient one."""
    cols = ["group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
            "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code", "auth_seq_id", "auth_comp_id",
            "auth_asym_id", "auth_atom_id", "B_iso_or_equiv", "occupancy", "label_entity_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num", "id"]
    hdr = ["data_pred", "#", "loop_"] + [f"_atom_site.{c} " for c in cols]
    one_to_three = {v: k for k, v in THREE_TO_ONE.items()}
    rows = []
    for i, (xyz, aa) in enumerate(zip(coords, seq), start=1):
        aa3 = one_to_three.get(aa, "UNK")
        rows.append(f"ATOM C CA . {aa3} A {i} . {i} {aa3} A CA 90.00 1.0 1 "
                    f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} 1 {i}")
    return "\n".join(hdr + rows + ["#", ""])


def _helix(nres, seed=0, noise=0.0):
    rng = np.random.default_rng(seed)
    t = np.arange(nres)
    c = np.stack([2.3 * np.cos(t * 1.75), 2.3 * np.sin(t * 1.75), 1.5 * t], 1)
    return c + (rng.normal(0, noise, c.shape) if noise else 0.0)


def self_test():
    import tempfile
    from types import SimpleNamespace
    print("SELF-TEST: fabricate two runs on disk and score them (no S3/creds/GPU).\n")
    ok = lambda m: print(f"  [ok] {m}")

    with tempfile.TemporaryDirectory() as td:
        NRES, NIDS = 60, 40
        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(NRES))
        base = _helix(NRES, seed=1)
        rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)   # 90 deg about z
        dira, dirb = os.path.join(td, "runA"), os.path.join(td, "runB")
        for d in (dira, dirb):
            os.makedirs(os.path.join(d, "structures"))
        ids = [str(1000 + i) for i in range(NIDS)]
        for n, oid in enumerate(ids):
            os.path.join(dira, "structures", f"orf{oid}.cif")
            with open(os.path.join(dira, "structures", f"orf{oid}.cif"), "w") as f:
                f.write(_fake_cif(base, seq))
            # B: rigid-body moved copy for most ids (should score TM ~ 1), progressively distorted
            # for the last few (should score lower) -- gives the aggregates a real spread to bin.
            cb = base @ rot.T + np.array([13.0, -4.0, 7.5])
            if n >= NIDS - 8:
                cb = cb + _helix(NRES, seed=2 + n, noise=1.0 + 0.9 * (n - (NIDS - 8)))- _helix(NRES, seed=2 + n)
            with open(os.path.join(dirb, "structures", f"orf{oid}.cif"), "w") as f:
                f.write(_fake_cif(cb, seq))
        # one id folded only in A -> must be reported, not crash
        with open(os.path.join(dira, "structures", "orf9999.cif"), "w") as f:
            f.write(_fake_cif(base, seq))

        # --- parsers agree ------------------------------------------------
        text = open(os.path.join(dira, "structures", f"orf{ids[0]}.cif")).read()
        cf, sf = parse_cif_fast(text)
        assert cf.shape == (NRES, 3) and sf == seq, (cf.shape, sf[:20])
        ok(f"fast parser: {NRES} CA + sequence recovered")
        try:
            cb2, sb2 = parse_cif_biopython(text)
            assert np.allclose(cf, cb2, atol=1e-3) and sf == sb2
            ok("biopython parser agrees with the fast parser")
        except ImportError:
            print("  [--] biopython not installed, skipped cross-check")

        # --- metric invariance --------------------------------------------
        r_id, tm_id = kabsch_rmsd_tm(base, base)
        assert r_id < 1e-9 and abs(tm_id - 1.0) < 1e-9, (r_id, tm_id)
        r_rt, tm_rt = kabsch_rmsd_tm(base, base @ rot.T + 5.0)
        assert r_rt < 1e-6 and abs(tm_rt - 1.0) < 1e-6, (r_rt, tm_rt)
        ok("kabsch TM: 1.0 on identical AND on rigid-body-moved copies (invariant)")

        def mkargs(**kw):
            d = dict(run_a=dira, run_b=dirb, label_a="base", label_b="msa",
                     out=os.path.join(td, "out"), ids=None, ids_from_listing=False,
                     no_intersect=False, cif_name="orf{id}.cif", parser="fast", metric="both",
                     workers=4, shard=None, limit=None, resume=True, sample=1000, seed=0,
                     flush_every=10 ** 9, s3_out=None)
            d.update(kw)
            return SimpleNamespace(**d)

        # --- full run -----------------------------------------------------
        a1 = mkargs()
        run(a1)
        rows = list(csv.DictReader(open(os.path.join(a1.out, "pairs.tsv")), delimiter="\t"))
        assert len(rows) == NIDS, f"expected {NIDS} rows, got {len(rows)}"
        assert all(r["seq_match"] == "exact" for r in rows), "same seq must compare as exact"
        ok(f"scored {len(rows)} pairs; worklist intersected the runs (orf9999 excluded: only in A)")

        rigid = [r for r in rows if float(r["tm_min"]) > 0.99]
        moved = [r for r in rows if float(r["tm_min"]) <= 0.99]
        assert len(rigid) == NIDS - 8, f"{len(rigid)} rigid-body pairs scored ~1 (want {NIDS-8})"
        assert len(moved) == 8, f"{len(moved)} distorted pairs scored <1 (want 8)"
        ok("TM-align: ~1.0 for inert pairs, <1 for distorted pairs")
        # kabsch vs tm-align, split by regime -- they are NOT related by a one-sided bound:
        #  * near-identical pairs (the regime --metric auto actually uses kabsch in): the two must
        #    agree tightly, but tm-align can land marginally BELOW kabsch, because its superposition
        #    is an iterative heuristic rather than a global optimum (observed -5e-3 on real lengths).
        #  * genuinely moved pairs: tm-align may find a materially BETTER correspondence than the
        #    fixed 1:1 one, so it legitimately exceeds kabsch (0.078 on these fixtures).
        gn = [float(r["tm_min"]) - float(r["tm_kabsch"]) for r in rows
              if r["tm_kabsch"] and float(r["tm_min"]) > 0.99]
        gf = [float(r["tm_min"]) - float(r["tm_kabsch"]) for r in rows
              if r["tm_kabsch"] and float(r["tm_min"]) <= 0.99]
        assert gn and max(abs(g) for g in gn) < 0.01, f"near-identical pairs must agree: {gn}"
        assert gf and min(gf) > -0.01, f"on moved pairs tm-align should not fall below kabsch: {min(gf)}"
        ok(f"kabsch vs tm-align: agree to {max(abs(g) for g in gn):.5f} on inert pairs; "
           f"tm-align up to {max(gf):+.4f} higher on moved pairs (better correspondence)")

        # --- aggregates are faithful + O(1) ------------------------------
        agg = json.load(open(os.path.join(a1.out, "aggregates.json")))
        assert agg["n_ok"] == NIDS == sum(agg["hist"]["tm"]), (agg["n_ok"], sum(agg["hist"]["tm"]))
        assert sum(agg["hist"]["rmsd"]) == NIDS and np.sum(agg["hist"]["tm_by_len"]) == NIDS
        tms = sorted(float(r["tm_min"]) for r in rows)
        assert abs(agg["moments"]["mean"] - float(np.mean(tms))) < 1e-9
        assert agg["tm_below"]["0.99"] == sum(1 for t in tms if t < 0.99)
        med = quantiles_from_hist(agg["hist"]["tm"], TM_BINS)["0.5"]
        assert abs(med - np.median(tms)) <= 0.005 + 1e-9, (med, np.median(tms))
        ok("aggregates: histograms sum to n, Welford mean exact, threshold counts exact, "
           "median within one bin")

        # --- missing structure is reported, not fatal ---------------------
        ids_file = os.path.join(td, "ids.txt")
        open(ids_file, "w").write("\n".join(ids[:3] + ["4242"]))
        a2 = mkargs(out=os.path.join(td, "out2"), ids=ids_file)
        run(a2)
        r2 = {r["orfid"]: r["status"] for r in
              csv.DictReader(open(os.path.join(a2.out, "pairs.tsv")), delimiter="\t")}
        assert r2["4242"] == "missing_a", r2
        assert sum(1 for v in r2.values() if v == "ok") == 3, r2
        ok("--ids honoured; an id absent from a run -> status=missing_a (no crash)")

        # --- 'auto' picks the right metric per pair -----------------------
        # Same sequence -> correspondence known -> kabsch (exact, ~1000x cheaper).
        # Truncated in one run -> no 1:1 correspondence -> must fall back to tm-align.
        dirc, dird = os.path.join(td, "runC"), os.path.join(td, "runD")
        for d in (dirc, dird):
            os.makedirs(os.path.join(d, "structures"))
        with open(os.path.join(dirc, "structures", "orf7001.cif"), "w") as f:
            f.write(_fake_cif(base, seq))                        # same length both runs
        with open(os.path.join(dird, "structures", "orf7001.cif"), "w") as f:
            f.write(_fake_cif(base @ rot.T + 2.0, seq))
        with open(os.path.join(dirc, "structures", "orf7002.cif"), "w") as f:
            f.write(_fake_cif(base, seq))                        # truncated in run D
        with open(os.path.join(dird, "structures", "orf7002.cif"), "w") as f:
            f.write(_fake_cif(base[:NRES - 12], seq[:NRES - 12]))
        a_auto = mkargs(run_a=dirc, run_b=dird, out=os.path.join(td, "auto"), metric="auto")
        run(a_auto)
        with open(os.path.join(a_auto.out, "pairs.tsv"), newline="") as fh:
            ra = {r["orfid"]: r for r in csv.DictReader(fh, delimiter="\t")}
        assert ra["7001"]["metric_used"] == "kabsch", ra["7001"]
        assert ra["7001"]["seq_match"] == "exact" and float(ra["7001"]["tm_min"]) > 0.999
        assert ra["7002"]["metric_used"] == "tmalign", ra["7002"]
        assert ra["7002"]["seq_match"] == "len_diff" and ra["7002"]["status"] == "ok"
        ok("--metric auto: kabsch when sequences match, tm-align fallback when one run truncated")

        # --- resume -------------------------------------------------------
        n_before = len(open(os.path.join(a2.out, "pairs.tsv")).readlines())
        run(mkargs(out=a2.out, ids=ids_file))
        n_after = len(open(os.path.join(a2.out, "pairs.tsv")).readlines())
        assert n_after == n_before, f"resume re-scored rows: {n_before} -> {n_after}"
        ok("resume: re-running appends nothing (ledger is the completion record)")

        # --- sharding: disjoint, complete, and merges back ----------------
        seen, N = [], 4
        for k in range(1, N + 1):
            ak = mkargs(out=os.path.join(td, "sh"), shard=f"{k}/{N}")
            run(ak)
            seen.append({r["orfid"] for r in csv.DictReader(
                open(os.path.join(ak.out, f"pairs.shard{k}of{N}.tsv")), delimiter="\t")})
        union = set().union(*seen)
        assert union == set(ids), f"shards missed {set(ids) - union} / extra {union - set(ids)}"
        for i in range(N):
            for j in range(i + 1, N):
                assert not (seen[i] & seen[j]), f"shards {i+1},{j+1} overlap: {seen[i] & seen[j]}"
        ok(f"--shard: {N} shards are disjoint and cover every pair exactly once")

        merge_mode(mkargs(out=os.path.join(td, "sh")))
        m = json.load(open(os.path.join(td, "sh", "aggregates.merged.json")))
        assert m["n_ok"] == agg["n_ok"], (m["n_ok"], agg["n_ok"])
        assert m["hist"]["tm"] == agg["hist"]["tm"], "merged histogram != single-run histogram"
        assert abs(m["moments"]["mean"] - agg["moments"]["mean"]) < 1e-9
        assert abs(m["moments"]["var"] - agg["moments"]["var"]) < 1e-9
        ok("--merge: 4 shard aggregates reassemble bin-for-bin into the single-run aggregate")

    print("\nSELF-TEST PASSED")
    return 0


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-a", help="reference run: s3://bucket/prefix/ or local dir (baseline)")
    ap.add_argument("--run-b", help="comparison run: s3://bucket/prefix/ or local dir (e.g. MSA)")
    ap.add_argument("--out", default="pairtm", help="output dir for ledger + aggregates (default pairtm/)")
    ap.add_argument("--label-a", default="", help="short name for run A in plots/summaries")
    ap.add_argument("--label-b", default="", help="short name for run B")
    ap.add_argument("--ids", help="file of ORFids to score (one per token; 'orf' prefix optional). "
                                  "Cheapest + deterministic; skips listing both runs.")
    ap.add_argument("--ids-from-listing", action="store_true",
                    help="build the work list by LISTING both runs' structures/ instead of reading "
                         "their tracker.csv (use if the trackers are absent)")
    ap.add_argument("--no-intersect", action="store_true",
                    help="with --ids, do not warn about ids absent from a run (they become missing_* rows)")
    ap.add_argument("--cif-name", default="orf{id}.cif", dest="cif_name",
                    help="structure filename template inside structures/ (default orf{id}.cif)")
    ap.add_argument("--metric", choices=["auto", "tmalign", "kabsch", "both"], default="auto",
                    help="auto (DEFAULT: kabsch when the two sequences match -- the correspondence is "
                         "then known and kabsch is exact and ~1000x cheaper -- else tmalign) | tmalign "
                         "(always, ~cubic in length) | kabsch (always; needs equal lengths) | both "
                         "(compute both and record the gap, for the one-off proxy justification)")
    ap.add_argument("--parser", choices=list(PARSERS), default="fast",
                    help="fast (default, ~1ms) | biopython (strict, ~20ms)")
    ap.add_argument("--workers", type=int, default=16,
                    help="thread pool size; the work is S3-latency-bound (default 16)")
    ap.add_argument("--shard", default=None, metavar="K/N",
                    help="score only band K of N (1-based), for fleet fan-out; falls back to "
                         "$PAIRTM_SHARD. Each shard writes its own pairs.shardKofN.tsv + "
                         "aggregates.shardKofN.json, so workers never clobber each other. "
                         "Combine them afterwards with --merge.")
    ap.add_argument("--limit", type=int, help="only score the first N pairs (smoke test)")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                    help="rescore ids already in the ledger (default: skip them)")
    ap.add_argument("--sample", type=int, default=20000,
                    help="reservoir subsample size kept in the aggregates for scatter/outliers")
    ap.add_argument("--seed", type=int, default=0, help="reservoir seed (reproducible subsample)")
    ap.add_argument("--flush-every", type=int, default=500,
                    help="flush ledger + rewrite aggregates every N pairs (progress/crash safety)")
    ap.add_argument("--s3-out", metavar="s3://bucket/prefix/",
                    help="also mirror the ledger + aggregates here at the end (fleet runs)")
    ap.add_argument("--merge", action="store_true",
                    help="do not score: merge --out's aggregates.shard*.json into aggregates.merged.json")
    ap.add_argument("--self-test", action="store_true", help="run offline correctness checks and exit")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(self_test())
    if a.merge:
        sys.exit(merge_mode(a))
    if not (a.run_a and a.run_b):
        ap.error("need --run-a and --run-b (or --self-test / --merge)")
    need = ["numpy"] + (["tmtools"] if a.metric in ("auto", "tmalign", "both") else [])
    need += ["boto3"] if str(a.run_a).startswith("s3://") or str(a.run_b).startswith("s3://") else []
    for m in need:
        try:
            __import__(m)
        except ImportError:
            sys.exit(f"missing dependency '{m}'.  pip install {' '.join(need)}")
    sys.exit(run(a))


if __name__ == "__main__":
    main()
