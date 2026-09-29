"""
molprobity_score.py -- figures Bb/Bc: MolProbity all-atom validation of predicted structures.

Runs ON the EC2 VM that has MolProbity built (see the notebook's VM setup section). Streams predicted
structures from S3, scores each one, and appends a row per structure to a resumable TSV ledger:

    molprobity.python molprobity_score.py \
        --run base=s3://petadex-protein-structures/esmfold2_paramsweep/s100_l20/ \
        --run msa=s3://petadex-protein-structures/esmfold2_paramsweep/s32_l5/ \
        --out mp_base_vs_msa --workers 16

`--run LABEL=SOURCE` is repeatable; SOURCE is an `s3://bucket/prefix/` or a local dir laid out as
esmfold2_local_predictor.py writes it (`<root>/structures/orf<ORFid>.cif`).

WHY THIS RUNS ON A VM AND NOT LOCALLY
  MolProbity's three components each need reference data that conda-forge's cctbx-base does NOT ship
  (CCP4 monomer library, chem_data/rotarama_data, chem_data/chemical_components), and clashscore also
  needs hydrogens added by `reduce`. The MolProbity bootstrap build pulls all of it at once, which is
  why the canonical score is produced on a purpose-built image rather than in the analysis env.

WHAT IS BEING MEASURED
  MolProbity score combines three terms; LOWER IS BETTER:

      MPscore = 0.426*ln(1+clashscore)
              + 0.330*ln(1+max(0, rota_outliers% - 1))
              + 0.250*ln(1+max(0, (100-rama_favored%) - 2))
              + 0.5

  * clashscore   = serious steric overlaps (>0.4 A) per 1000 atoms. NEEDS HYDROGENS -- our predictions
                   have none, so `reduce` adds them first. This is the term that usually dominates for
                   predicted structures.
  * rota/rama    = sidechain rotamer outliers and backbone Ramachandran favoured fraction.
  All three are required: there is no partial MolProbity score.

  These are ABSOLUTE geometry measures, not accuracy. A good MolProbity score means the model is
  physically plausible, NOT that it matches the real structure.

ENGINES (--engine)
  api  in-process cctbx call. Preferred: cctbx takes 1-3 s just to START, which would otherwise
       dominate a 1-3 s job, so per-structure subprocesses are several times slower.
  cli  shells out to `molprobity.molprobity` per structure and parses its JSON. Slower, but depends
       only on the documented CLI surface.
  auto (default) try api once, fall back to cli.

  RUN `--check-env` FIRST ON THE VM. It reports which engine works and scores one structure, so the
  engine question is settled by measurement before a batch is launched.

Deps on the VM: a MolProbity build (provides molprobity.python + molprobity.molprobity) + boto3.
"""
import argparse
import csv
import io
import math
import os
import sys
import tempfile
import time
import threading
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

# ---------------------------------------------------------------------------
# the score
# ---------------------------------------------------------------------------
def molprobity_score(clashscore, rota_outliers, rama_favored):
    """Published MolProbity score (Williams et al. 2018). Lower is better; 0.5 is the floor for a
    structure with no clashes, no rotamer outliers and 100% favoured Ramachandran."""
    if None in (clashscore, rota_outliers, rama_favored):
        return None
    rama_iffy = 100.0 - rama_favored
    return max(0.0,
               0.426 * math.log(1 + clashscore)
               + 0.330 * math.log(1 + max(0.0, rota_outliers - 1.0))
               + 0.250 * math.log(1 + max(0.0, rama_iffy - 2.0))
               + 0.5)


# ---------------------------------------------------------------------------
# structure source: S3 stream (default) or local dir
# ---------------------------------------------------------------------------
class StructSource:
    """Reads <root>/structures/<cif_name>. `root` is s3://bucket/prefix/ or a local dir.

    One boto3 client PER THREAD -- clients are only nominally thread-safe and sharing one across a
    pool serialises on its connection pool.
    """

    def __init__(self, label, root, cif_name="orf{id}.cif", retries=4):
        self.label, self.root, self.cif_name, self.retries = label, root, cif_name, retries
        self.is_s3 = str(root).startswith("s3://")
        self._tl = threading.local()
        if self.is_s3:
            self.bucket, _, self.prefix = root[len("s3://"):].partition("/")
            self.prefix = self.prefix.strip("/")

    def __getstate__(self):
        # threading.local() cannot be pickled, and a boto3 client must never be inherited across a
        # fork (shared sockets/SSL state corrupt). Drop it; the property rebuilds it on demand.
        d = self.__dict__.copy()
        d.pop("_tl", None)
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._tl = threading.local()

    @property
    def s3(self):
        if getattr(self._tl, "c", None) is None:
            import boto3
            self._tl.c = boto3.client("s3")
        return self._tl.c

    def get_text(self, oid):
        """CIF text, or None if this run has no structure for `oid`."""
        name = self.cif_name.format(id=oid)
        if not self.is_s3:
            for p in (os.path.join(self.root, "structures", name), os.path.join(self.root, name)):
                if os.path.exists(p):
                    with open(p) as f:
                        return f.read()
            return None
        key = "/".join(p for p in (self.prefix, "structures", name) if p)
        for attempt in range(self.retries):
            try:
                return self.s3.get_object(Bucket=self.bucket, Key=key)["Body"].read().decode()
            except Exception as e:                        # noqa: BLE001
                code = str(getattr(e, "response", {}).get("Error", {}).get("Code", ""))
                if code in ("NoSuchKey", "404", "NoSuchBucket") or type(e).__name__ == "NoSuchKey":
                    return None
                if attempt == self.retries - 1:
                    raise
                time.sleep(0.4 * 2 ** attempt)
        return None

    def list_ids(self):
        """Every ORFid with a structure under this root."""
        import re
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
            kw = {"Bucket": self.bucket,
                  "Prefix": "/".join(p for p in (self.prefix, "structures/") if p)}
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


def norm_id(x):
    """Canonical bare ORFid ('orf3772973' and '3772973' must land on the same key)."""
    x = str(x).strip()
    return x[3:] if x.startswith("orf") else x


# ---------------------------------------------------------------------------
# engines
# ---------------------------------------------------------------------------
def _get(obj, name):
    """cctbx exposes some of these as methods and some as attributes, and it has moved between the
    two across versions. Try both rather than pin one and break on upgrade."""
    v = getattr(obj, name, None)
    if v is None:
        return None
    try:
        return v() if callable(v) else v
    except Exception:                                     # noqa: BLE001
        return None


def score_api(path):
    """In-process cctbx validation. Returns a dict of the three components + n_res."""
    from iotbx.data_manager import DataManager
    from mmtbx.validation.molprobity import molprobity as molprobity_cls

    dm = DataManager()
    dm.process_model_file(path)
    model = dm.get_model()
    try:
        model.process(make_restraints=True)               # reduce/probe need restraints
    except TypeError:
        model.process()                                   # older signature
    mp = molprobity_cls(model=model)

    n_res = None
    try:
        n_res = len(list(model.get_hierarchy().residue_groups()))
    except Exception:                                     # noqa: BLE001
        pass
    return {
        "clashscore": _num(_get(mp, "clashscore")),
        "pct_rota_outliers": _num(_get(mp, "rota_outliers")),
        "pct_rama_favored": _num(_get(mp, "rama_favored")),
        "pct_rama_outliers": _num(_get(mp, "rama_outliers")),
        "molprobity_score_native": _num(_get(mp, "molprobity_score")),
        "n_res": n_res,
    }


def score_cli(path, timeout=600):
    """Shell out to `molprobity.molprobity` and parse its JSON."""
    import json
    import subprocess
    cmd = ["molprobity.molprobity", path, "output.quiet=True", "json=True"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"molprobity.molprobity exited {r.returncode}: {r.stderr[-400:]}")
    txt = r.stdout.strip()
    start = txt.find("{")
    if start < 0:
        raise RuntimeError(f"no JSON in molprobity output: {txt[-400:]}")
    d = json.loads(txt[start:])
    flat = _flatten(d)

    def pick(*names):
        for n in names:
            if n in flat:
                return _num(flat[n])
        return None

    return {
        "clashscore": pick("clashscore", "clash.clashscore", "summary_results.clashscore"),
        "pct_rota_outliers": pick("rota_outliers", "rotamer_outliers", "rota.outliers"),
        "pct_rama_favored": pick("rama_favored", "ramachandran_favored", "rama.favored"),
        "pct_rama_outliers": pick("rama_outliers", "ramachandran_outliers", "rama.outliers"),
        "molprobity_score_native": pick("molprobity_score", "MolProbity_score"),
        "n_res": pick("n_residues", "num_residues"),
    }


def _flatten(d, prefix=""):
    """Flatten nested JSON to dotted keys, keeping the leaf name as an alias too -- the CLI's exact
    nesting varies by version, so match on leaf names rather than a hardcoded path."""
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            if isinstance(v, (dict, list)):
                out.update(_flatten(v, key))
            else:
                out[key] = v
                out.setdefault(str(k), v)
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(_flatten(v, f"{prefix}.{i}"))
    return out


def _num(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


ENGINES = {"api": score_api, "cli": score_cli}


# ---------------------------------------------------------------------------
# ledger
# ---------------------------------------------------------------------------
COLS = ["orfid", "run", "n_res", "clashscore", "pct_rota_outliers", "pct_rama_favored",
        "pct_rama_outliers", "molprobity_score", "molprobity_score_native", "engine",
        "wall_s", "status", "error"]


def read_done(path):
    """(run, orfid) pairs already scored, so a restart resumes instead of redoing work."""
    if not os.path.exists(path):
        return set()
    done = set()
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("orfid") and row.get("status") == "ok":
                done.add((row.get("run", ""), row["orfid"]))
    return done


# ---------------------------------------------------------------------------
# scoring one structure
# ---------------------------------------------------------------------------
def score_one(oid, src, engine_name, engine_fn):
    rec = {c: "" for c in COLS}
    rec.update(orfid=oid, run=src.label, engine=engine_name, status="error")
    t0 = time.perf_counter()

    text = src.get_text(oid)
    if text is None:
        rec["status"] = "missing"
        return rec

    # cctbx and the CLI both want a path; keep the temp file out of the way and always clean up.
    fd, tmp = tempfile.mkstemp(suffix=".cif")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
        d = engine_fn(tmp)
    except Exception as e:                                # noqa: BLE001  one bad structure must not kill the run
        rec["status"] = "score_fail"
        rec["error"] = f"{type(e).__name__}: {e}"[:300]
        rec["wall_s"] = round(time.perf_counter() - t0, 2)
        return rec
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass

    mp = molprobity_score(d.get("clashscore"), d.get("pct_rota_outliers"), d.get("pct_rama_favored"))
    rec.update({k: ("" if d.get(k) is None else d[k]) for k in
                ("n_res", "clashscore", "pct_rota_outliers", "pct_rama_favored",
                 "pct_rama_outliers")})
    rec["molprobity_score_native"] = "" if d.get("molprobity_score_native") is None \
        else round(d["molprobity_score_native"], 3)
    rec["molprobity_score"] = "" if mp is None else round(mp, 3)
    rec["wall_s"] = round(time.perf_counter() - t0, 2)
    rec["status"] = "ok" if mp is not None else "incomplete"
    if mp is None:
        rec["error"] = "a component was missing -- no partial MolProbity score is defined"
    return rec


# ---------------------------------------------------------------------------
# sharding (same 'K/N' contract as the predictor and pair_tm.py)
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


def select_shard(items, k, n):
    """Contiguous band K of N over a deterministic ordering -> shards never overlap or gap."""
    ordered = sorted(items, key=lambda s: (len(str(s)), str(s)))
    b = [round(i * len(ordered) / n) for i in range(n + 1)]
    return ordered[b[k - 1]:b[k]]


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def pick_engine(name, sample_path=None):
    """Resolve --engine. 'auto' prefers the in-process API and falls back to the CLI."""
    if name in ENGINES:
        # Warm the same imports the 'auto' path below does. An explicit --engine api would
        # otherwise leave the first cctbx import to happen inside N worker threads at once,
        # which deadlocks on the import lock (this bit pair_tm.py in 6a).
        if name == "api":
            from iotbx.data_manager import DataManager           # noqa: F401
            from mmtbx.validation.molprobity import molprobity    # noqa: F401
        if sample_path:
            ENGINES[name](sample_path)
        return name, ENGINES[name]
    order = ["api", "cli"]
    errs = []
    for n in order:
        try:
            if n == "api":
                from iotbx.data_manager import DataManager           # noqa: F401
                from mmtbx.validation.molprobity import molprobity   # noqa: F401
            else:
                import subprocess
                subprocess.run(["molprobity.molprobity", "--help"],
                               capture_output=True, timeout=120)
            if sample_path:
                ENGINES[n](sample_path)
            return n, ENGINES[n]
        except Exception as e:                            # noqa: BLE001
            errs.append(f"  {n}: {type(e).__name__}: {e}")
    raise SystemExit("no working MolProbity engine:\n" + "\n".join(errs)
                     + "\n(are you running under `molprobity.python` on the built VM?)")


def build_worklist(a, sources):
    if a.ids:
        with open(a.ids) as f:
            return sorted({norm_id(x) for x in f.read().split() if x.strip()}), f"--ids {a.ids}"
    sets = [{norm_id(x) for x in s.list_ids()} for s in sources]
    if a.intersect and len(sets) > 1:
        ids = set.intersection(*sets)
        extra = sum(len(s - ids) for s in sets)
        return sorted(ids), f"listing, intersected across {len(sources)} runs ({extra} unpaired dropped)"
    return sorted(set.union(*sets)), "listing (union)"


# --- worker-process state --------------------------------------------------
# cctbx is Boost.Python, which holds the GIL through its C++ work, so a thread pool serialises
# almost completely (measured: 16 threads == 4 threads == ~4.5 s/structure). Processes get real
# parallelism. Everything a worker needs is rebuilt inside the worker rather than pickled across,
# which keeps boto3 clients per-process and avoids inheriting sockets over a fork.
_W = {}


def _init_worker(engine_name, specs, cif_name):
    """Runs once per worker process (or once in-process for --pool thread)."""
    global _W
    name, fn = pick_engine(engine_name)          # imports cctbx HERE, before any scoring
    _W = {"engine_name": name, "engine_fn": fn,
          "sources": {label: StructSource(label, root, cif_name) for label, root in specs}}


def _work(job):
    oid, label = job
    try:
        return score_one(oid, _W["sources"][label], _W["engine_name"], _W["engine_fn"])
    except Exception as e:                                # noqa: BLE001
        r = {c: "" for c in COLS}
        r.update(orfid=oid, run=label, engine=_W.get("engine_name", ""), status="error",
                 error=f"{type(e).__name__}: {e}"[:300])
        return r


def fmt_hms(seconds):
    s = int(max(seconds, 0))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f"{h}:{m:02d}:{sec:02d}" if h else f"{m}:{sec:02d}"


class Progress:
    """Dependency-free progress bar (no tqdm on the MolProbity conda env).

    Redraws in place on a tty; falls back to one line every `every` structures when stdout is
    redirected, so `> run.log` and `nohup` stay readable instead of accumulating \\r spam.
    """

    def __init__(self, total, width=32, every=25, stream=None):
        self.total = max(int(total), 1)
        self.width = width
        self.every = max(int(every), 1)
        self.stream = stream if stream is not None else sys.stdout
        self.tty = hasattr(self.stream, "isatty") and self.stream.isatty()
        self.t0 = time.perf_counter()

    def _stats(self, n):
        elapsed = time.perf_counter() - self.t0
        rate = n / max(elapsed, 1e-9)
        eta = (self.total - n) / rate if rate > 0 else 0.0
        return elapsed, rate, eta

    def update(self, n, ok, bad=0):
        elapsed, rate, eta = self._stats(n)
        pct = 100.0 * n / self.total
        tail = (f"{n}/{self.total} {pct:5.1f}%  ok={ok}" + (f" bad={bad}" if bad else "")
                + f"  {rate:.2f}/s  {fmt_hms(elapsed)}<{fmt_hms(eta)}")
        if self.tty:
            filled = int(self.width * n / self.total)
            bar = "#" * filled + "-" * (self.width - filled)
            self.stream.write("\r  [" + bar + "] " + tail.ljust(58)[:58])
            self.stream.flush()
        elif n % self.every == 0 or n == self.total:
            self.stream.write("  " + tail + "\n")
            self.stream.flush()

    def close(self):
        if self.tty:
            self.stream.write("\n")
            self.stream.flush()


def run(a):
    shard = parse_shard(a.shard if a.shard is not None else os.environ.get("MOLPROBITY_SHARD"))
    tag = f".shard{shard[0]}of{shard[1]}" if shard else ""
    os.makedirs(a.out, exist_ok=True)
    ledger = os.path.join(a.out, f"molprobity{tag}.tsv")

    sources = []
    for spec in a.run:
        if "=" not in spec:
            raise SystemExit(f"--run must be LABEL=SOURCE, got {spec!r}")
        label, _, root = spec.partition("=")
        sources.append(StructSource(label.strip(), root.strip(), a.cif_name))
    for s in sources:
        print(f"run '{s.label}'  {s.root}")

    ids, how = build_worklist(a, sources)
    print(f"{len(ids)} ORFs (work list from {how})")
    if shard:
        before = len(ids)
        ids = select_shard(ids, *shard)
        print(f"shard {shard[0]}/{shard[1]}: {len(ids)}/{before} ORFs")
    if a.limit:
        ids = ids[:a.limit]
        print(f"--limit {a.limit}: {len(ids)} ORFs")

    jobs = [(oid, s.label) for s in sources for oid in ids]
    done = read_done(ledger) if a.resume else set()
    if done:
        jobs = [(o, lab) for (o, lab) in jobs if (lab, o) not in done]
        print(f"resume: {len(done)} already scored, {len(jobs)} left")
    if not jobs:
        print("nothing to do -- ledger already complete for this shard")
        return 0

    # Settle the engine ONCE, in the main process, before any workers start: it imports cctbx (or
    # probes the CLI), and racing that from N threads is exactly the kind of first-import contention
    # that deadlocks.
    engine_name, engine_fn = pick_engine(a.engine)
    print(f"engine: {engine_name}   pool: {a.pool}   workers: {a.workers}\n")

    initargs = (engine_name, [(s.label, s.root) for s in sources], a.cif_name)
    if a.pool == "process":
        make_pool = lambda: ProcessPoolExecutor(                       # noqa: E731
            max_workers=a.workers, initializer=_init_worker, initargs=initargs)
    else:
        _init_worker(*initargs)          # threads share one interpreter: initialise once, here
        make_pool = lambda: ThreadPoolExecutor(max_workers=a.workers)  # noqa: E731

    lock = threading.Lock()
    counts = Counter()
    t0 = time.perf_counter()
    fresh = not os.path.exists(ledger) or os.path.getsize(ledger) == 0
    with open(ledger, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        if fresh:
            w.writeheader()

        n = 0
        prog = Progress(len(jobs), every=a.flush_every)
        with make_pool() as pool:
            # chunksize=1 keeps results (and the bar) flowing one structure at a time rather than
            # arriving in bursts -- at ~4 s each the scheduling overhead is irrelevant.
            for rec in pool.map(_work, jobs, chunksize=1):
                with lock:
                    w.writerow(rec)
                    counts[rec["status"]] += 1
                    n += 1
                    if n % a.flush_every == 0:       # durability is on its own clock now,
                        fh.flush()                   # decoupled from how often we redraw
                    prog.update(n, counts["ok"], n - counts["ok"])
        prog.close()

    elapsed = time.perf_counter() - t0
    print(f"\nscored {n} structures in {elapsed:.1f}s ({n/max(elapsed,1e-9):.2f}/s)")
    print("  " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"  ledger -> {ledger}")
    if counts.get("ok"):
        summarise(ledger)
    if a.s3_out:
        mirror(a.s3_out, [ledger])
    return 0


def summarise(ledger):
    """Per-run medians -- a sanity read before pulling the TSV down for the figures."""
    import statistics
    by = {}
    with open(ledger, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["status"] != "ok":
                continue
            v = _num(row["molprobity_score"])
            c = _num(row["clashscore"])
            if v is not None:
                by.setdefault(row["run"], []).append((v, c))
    print()
    for label, vals in sorted(by.items()):
        mp = [x[0] for x in vals]
        cs = [x[1] for x in vals if x[1] is not None]
        print(f"  {label:<10} n={len(mp):<5} MolProbity median={statistics.median(mp):.3f} "
              f"mean={statistics.mean(mp):.3f}"
              + (f"   clashscore median={statistics.median(cs):.2f}" if cs else ""))
    print("  (lower is better; these are geometry, NOT accuracy)")


def mirror(s3_out, paths):
    import boto3
    c = boto3.client("s3")
    bucket, _, prefix = s3_out[len("s3://"):].partition("/")
    for p in paths:
        key = "/".join(x for x in (prefix.strip("/"), os.path.basename(p)) if x)
        with open(p, "rb") as f:
            c.put_object(Bucket=bucket, Key=key, Body=f.read())
        print(f"  mirrored -> s3://{bucket}/{key}")


# ---------------------------------------------------------------------------
# --check-env : settle the engine question by measurement, before a batch
# ---------------------------------------------------------------------------
def check_env(a):
    print("MolProbity environment check")
    print("=" * 60)
    print(f"python: {sys.version.split()[0]}  ({sys.executable})")
    for mod in ("iotbx", "mmtbx", "boto3"):
        try:
            m = __import__(mod)
            print(f"  OK    import {mod}  ({getattr(m, '__file__', '?')})")
        except Exception as e:                            # noqa: BLE001
            print(f"  FAIL  import {mod}: {type(e).__name__}: {e}")

    path = a.check_file
    if not path:
        if not a.run:
            print("\n(pass --check-file <local.cif> or --run LABEL=SOURCE to score a real structure)")
            return 0
        label, _, root = a.run[0].partition("=")
        src = StructSource(label, root, a.cif_name)
        ids = sorted(src.list_ids())
        if not ids:
            print(f"\nno structures found under {root}")
            return 1
        text = src.get_text(ids[0])
        fd, path = tempfile.mkstemp(suffix=".cif")
        with os.fdopen(fd, "w") as f:
            f.write(text)
        print(f"\nfetched {ids[0]} from {root}")

    print("\nengine trials on one structure:")
    ok_any = False
    for name in ("api", "cli"):
        t = time.perf_counter()
        try:
            d = ENGINES[name](path)
            mp = molprobity_score(d.get("clashscore"), d.get("pct_rota_outliers"),
                                  d.get("pct_rama_favored"))
            print(f"  OK    {name:<4} {time.perf_counter()-t:6.2f}s  clash={d.get('clashscore')} "
                  f"rota_out={d.get('pct_rota_outliers')} rama_fav={d.get('pct_rama_favored')} "
                  f"-> MPscore={None if mp is None else round(mp,3)} "
                  f"(native {d.get('molprobity_score_native')})")
            ok_any = True
        except Exception as e:                            # noqa: BLE001
            print(f"  FAIL  {name:<4} {time.perf_counter()-t:6.2f}s  {type(e).__name__}: {e}")
            if a.verbose:
                traceback.print_exc()
    print("\n" + ("pick the faster working engine with --engine" if ok_any else
                  "NO ENGINE WORKS -- is this running under `molprobity.python` on the built VM?"))
    return 0 if ok_any else 1


# ---------------------------------------------------------------------------
# self-test -- everything EXCEPT the cctbx call (which needs the VM)
# ---------------------------------------------------------------------------
def self_test():
    import shutil
    print("SELF-TEST: formula, ledger/resume, sharding, source, TSV. "
          "(The cctbx scoring call needs the VM -- use --check-env there.)\n")
    ok = lambda m: print(f"  [ok] {m}")

    # --- formula ----------------------------------------------------------
    assert abs(molprobity_score(0.0, 0.0, 100.0) - 0.5) < 1e-12, "perfect structure must score 0.5"
    assert molprobity_score(0.0, 1.0, 98.0) == molprobity_score(0.0, 0.0, 100.0), \
        "the 1% rotamer and 2% rama allowances must be free"
    worse = [molprobity_score(c, 0, 100) for c in (0, 1, 5, 20, 50)]
    assert all(x < y for x, y in zip(worse, worse[1:])), f"must increase with clashscore: {worse}"
    assert molprobity_score(5, 10, 90) > molprobity_score(5, 0, 100)
    assert molprobity_score(None, 0, 100) is None, "a missing component must not fake a score"
    ok("formula: 0.5 floor, allowances free, monotonic in clashscore, None on missing component")

    # --- source + ledger + resume + shard --------------------------------
    td = tempfile.mkdtemp()
    try:
        runs = {}
        for label in ("base", "msa"):
            d = os.path.join(td, label)
            os.makedirs(os.path.join(d, "structures"))
            for i in range(12):
                with open(os.path.join(d, "structures", f"orf{3700000+i}.cif"), "w") as f:
                    f.write("data_x\n")
            runs[label] = d
        # one ORF only in 'base' -> intersect must drop it, union must keep it
        with open(os.path.join(runs["base"], "structures", "orf9999999.cif"), "w") as f:
            f.write("data_x\n")

        sa = StructSource("base", runs["base"])
        sb = StructSource("msa", runs["msa"])
        assert len(sa.list_ids()) == 13 and len(sb.list_ids()) == 12
        assert sa.get_text("3700000") is not None and sa.get_text("404404") is None
        ok("source: lists ids, reads a structure, returns None for a missing one")

        class A:
            ids = None
            intersect = True
        inter, _ = build_worklist(A, [sa, sb])
        A.intersect = False
        union, _ = build_worklist(A, [sa, sb])
        assert len(inter) == 12 and "9999999" not in inter, inter
        assert len(union) == 13 and "9999999" in union
        ok("worklist: --intersect drops unpaired ORFs (Bc needs pairs); union keeps them")

        led = os.path.join(td, "led.tsv")
        with open(led, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
            w.writeheader()
            w.writerow({"orfid": "1", "run": "base", "status": "ok", "molprobity_score": 1.5})
            w.writerow({"orfid": "2", "run": "base", "status": "score_fail"})
        done = read_done(led)
        assert ("base", "1") in done, done
        assert ("base", "2") not in done, "a failed row must be RETRIED, not treated as done"
        ok("resume: ok rows are skipped, failed rows are retried")

        ids = [str(i) for i in range(97)]
        seen, N = [], 5
        for k in range(1, N + 1):
            seen.append(set(select_shard(ids, k, N)))
        assert set().union(*seen) == set(ids), "shards must cover every id"
        for i in range(N):
            for j in range(i + 1, N):
                assert not (seen[i] & seen[j]), f"shards {i+1},{j+1} overlap"
        ok(f"--shard: {N} shards are disjoint and cover all {len(ids)} ids exactly once")

        # --- a failing engine must be recorded, not fatal ------------------
        def boom(_):
            raise RuntimeError("simulated cctbx failure")
        rec = score_one("3700000", sa, "api", boom)
        assert rec["status"] == "score_fail" and "simulated" in rec["error"], rec
        assert rec["molprobity_score"] == "", "a failed structure must not carry a score"
        ok("a structure that fails to score -> status=score_fail, no score, run continues")

        rec = score_one("404404", sa, "api", boom)
        assert rec["status"] == "missing", rec
        ok("a structure absent from the run -> status=missing (no crash)")

        # --- a partial result must not produce a score ---------------------
        rec = score_one("3700000", sa, "api", lambda p: {"clashscore": 2.0,
                                                         "pct_rota_outliers": None,
                                                         "pct_rama_favored": 96.0})
        assert rec["status"] == "incomplete" and rec["molprobity_score"] == "", rec
        ok("missing component -> status=incomplete, NOT a partial MolProbity score")

        # --- CLI JSON shape tolerance --------------------------------------
        flat = _flatten({"summary_results": {"a.pdb": {"clashscore": 3.5, "rota_outliers": 1.2}}})
        assert flat["clashscore"] == 3.5 and flat["rota_outliers"] == 1.2, flat
        ok("CLI JSON: leaf-name matching survives the report's nesting")
    finally:
        shutil.rmtree(td, ignore_errors=True)

    print("\nSELF-TEST PASSED")
    print("NOTE: the cctbx scoring call itself is NOT covered here -- run --check-env on the VM.")
    return 0


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", default=[], metavar="LABEL=SOURCE",
                    help="repeatable: a run to score, e.g. base=s3://bucket/prefix/ (or a local dir)")
    ap.add_argument("--out", default="molprobity", help="output dir for the ledger (default molprobity/)")
    ap.add_argument("--engine", default="auto", choices=["auto", "api", "cli"],
                    help="auto (default: in-process api, else cli) | api | cli")
    ap.add_argument("--ids", help="file of ORFids to score (one per token; 'orf' prefix optional)")
    ap.add_argument("--no-intersect", dest="intersect", action="store_false",
                    help="score every ORF in each run instead of only those present in ALL runs "
                         "(the paired Bc figure wants the intersection, which is the default)")
    ap.add_argument("--cif-name", default="orf{id}.cif", dest="cif_name",
                    help="structure filename template inside structures/ (default orf{id}.cif)")
    ap.add_argument("--workers", type=int, default=8,
                    help="pool size; scoring is CPU-bound, so match this to vCPUs (default 8)")
    ap.add_argument("--pool", default="process", choices=["process", "thread"],
                    help="process (default) gives real parallelism -- cctbx is Boost.Python and "
                         "holds the GIL, so threads serialise. 'thread' is a debugging fallback.")
    ap.add_argument("--shard", default=None, metavar="K/N",
                    help="score only band K of N (1-based); falls back to $MOLPROBITY_SHARD. Each "
                         "shard writes its own molprobity.shardKofN.tsv")
    ap.add_argument("--limit", type=int, help="only score the first N ORFs (smoke test)")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                    help="rescore ORFs already in the ledger (default: skip successful ones)")
    ap.add_argument("--flush-every", type=int, default=25,
                    help="flush the ledger + print progress every N structures")
    ap.add_argument("--s3-out", metavar="s3://bucket/prefix/", help="mirror the ledger here at the end")
    ap.add_argument("--check-env", action="store_true",
                    help="RUN THIS FIRST ON THE VM: report which engine works and score one structure")
    ap.add_argument("--check-file", help="local .cif for --check-env (else the first ORF of --run)")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--self-test", action="store_true",
                    help="offline checks of everything except the cctbx call")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(self_test())
    if a.check_env:
        sys.exit(check_env(a))
    if not a.run:
        ap.error("need at least one --run LABEL=SOURCE (or --self-test / --check-env)")
    sys.exit(run(a))


if __name__ == "__main__":
    main()
