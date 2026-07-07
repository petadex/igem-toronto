"""
Fold every 60pid family centroid with ESMFold2 and benchmark the run.

Model:  esmfold2-fast-2026-05  (https://biohub.ai/models/esmfold2)
        Folding uses the biohub-hosted ESMFold2 SDK: esm.sdk.esmfold2_client +
        FoldingConfig, input_builder.ProteinInput/StructurePredictionInput, and
        client.fold_all_atom(...). Inference is REMOTE (runs on biohub's GPUs), so
        this box needs no GPU -- it is an API client + S3 writer.

Input:  a CSV whose rows are the centroid amino-acid sequences to fold (produced by
        the SQL/Athena step, schema TBD). Columns are auto-detected; override with
        --seq-col / --id-col. If no id column exists, ids are synthesised as row<N>.

Output (mirrors resources/260530_issue79_esmfold_db_download so structures drop
        straight into the existing Atlas-derived pipeline):
          <outdir>/structures/<id>.cif   all-atom mmCIF (SDK .complex.to_mmcif())
          <outdir>/arrays/<id>.npz       per_residue_plddt[L], pae[L,L], residue_index[L]
          <outdir>/metrics.csv           one row per (centroid, arm): timing + confidence
        structures/, arrays/ and metrics.csv are uploaded to the default S3 sink
        (s3://petadex-protein-structures/esmfold2-centroids/; --s3 '' to disable).
        metrics.csv doubles as the resume ledger (skips any (centroid_id, arm) already
        status=ok); it is mirrored to S3 every ~5 min + at exit and seeded back from S3
        on a fresh box, so resume + the benchmark table survive an ephemeral VM.

Benchmark arms:
    single : plain single-sequence folding (ESMFold's native mode).            [implemented]
    msa    : fold with the "petadex MSA+ESM" conditioning.                     [TODO — stub]
             The MSA source is not decided yet; build_msa() is a documented seam
             to be filled in a later session once the actual MSAs exist. Selecting
             --arm msa today raises NotImplementedError on purpose.

Install:  pip install esm@git+https://github.com/Biohub/esm.git@main numpy boto3 tqdm

Auth:  export BIOHUB_TOKEN=<your biohub.ai API token>   (or --token)

SDK calls are aligned to the official biohub ESMFold2 tutorial (fold_all_atom + the
result's .plddt / .ptm / .pae / .complex.to_mmcif()). The tutorial documents remote
inference only -- there is no local-GPU loader -- so --backend local is disabled.

Example:
    # smoke-test 5 centroids against the remote API, kept local (no upload)
    python esmfold2_predictor.py centroids.csv --limit 5 --s3 ''

    # full run; uploads to the default S3 sink
    python esmfold2_predictor.py centroids.csv --outdir out

    # faster/cheaper pass with fewer refinement loops
    python esmfold2_predictor.py centroids.csv --num-loops 10
"""
import os
import csv
import sys
import time
import hashlib
import argparse
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from itertools import islice

import numpy as np

try:
    from tqdm import tqdm
except ImportError:                                   # graceful no-op if tqdm is absent
    class tqdm:                                        # noqa: N801
        def __init__(self, iterable=None, **kw): self.iterable = iterable
        def __iter__(self): return iter(self.iterable or [])
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def update(self, n=1): pass
        def set_postfix(self, **kw): pass
        @staticmethod
        def write(msg, *a, **k): print(msg)

DEFAULT_MODEL = "esmfold2-fast-2026-05"
DEFAULT_URL   = "https://biohub.ai"
S3_DEST       = "s3://petadex-protein-structures/esmfold2-centroids/"   # canonical sink (issue #6)

# Column-name candidates for auto-detecting the CSV schema (case-insensitive).
SEQ_COLS = ("sequence", "seq", "aa_sequence", "protein_sequence", "aaseq", "aa", "prot_seq")
ID_COLS  = ("centroid_id", "centroid", "cluster_id", "cluster", "id", "name",
            "header", "accession", "rep", "representative")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")           # 20 + ambiguity/seleno codes


# --- fold result container -------------------------------------------------
@dataclass
class FoldResult:
    sequence: str
    cif_text: str                       # mmCIF from the SDK's .complex.to_mmcif()
    per_residue_plddt: np.ndarray       # [L] float (0-100)
    ptm: float
    pae: np.ndarray | None              # [L, L] float (Ångströms), or None if not requested
    iptm: float | None = None           # interface pTM (complexes only; None for single chain)

    @property
    def mean_plddt(self) -> float:
        return float(np.mean(self.per_residue_plddt)) if self.per_residue_plddt.size else float("nan")


def _np(x):
    """Coerce a torch tensor / list / scalar to a numpy array (detaching if needed)."""
    if x is None:
        return None
    if hasattr(x, "detach"):                          # torch.Tensor
        x = x.detach().cpu().numpy()
    return np.asarray(x)


def _extract_result(sequence, result) -> FoldResult:
    """Normalise a client.fold_all_atom(...) result into a FoldResult.

    Aligned to the official biohub ESMFold2 tutorial: the result exposes
    `.complex.to_mmcif()` (structure), `.plddt` (per-residue tensor, 0-1),
    `.ptm` (scalar), `.pae` ([L,L] tensor, present only when FoldingConfig
    include_pae=True), and `.iptm` (interface pTM, complexes only).
    """
    if isinstance(result, Exception):                 # SDK may hand back an error object
        raise result

    cif = result.complex.to_mmcif()

    plddt = _np(result.plddt).astype(np.float32).reshape(-1)
    if plddt.size and float(np.nanmax(plddt)) <= 1.0:   # 0-1 -> 0-100 (Atlas convention)
        plddt = plddt * 100.0

    ptm = float(_np(result.ptm))

    pae = getattr(result, "pae", None)                # absent unless include_pae=True
    if pae is not None:
        pae = _np(pae).astype(np.float32)
        if pae.ndim == 3 and pae.shape[0] == 1:
            pae = pae[0]

    iptm_v = getattr(result, "iptm", None)            # only meaningful for complexes
    iptm = float(_np(iptm_v)) if iptm_v is not None else None

    return FoldResult(sequence, cif, plddt, ptm, pae, iptm)


# --- backends --------------------------------------------------------------
class ApiBackend:
    """Remote inference via the biohub-hosted ESMFold2 client (esm.sdk.esmfold2_client)."""

    def __init__(self, model=DEFAULT_MODEL, url=DEFAULT_URL, token=None, config_kwargs=None):
        try:
            from esm.sdk import esmfold2_client
            from esm.sdk.api import FoldingConfig
        except ImportError as e:
            raise SystemExit("API backend needs the esm SDK: "
                             "pip install esm@git+https://github.com/Biohub/esm.git@main "
                             "(see https://biohub.ai/models/esmfold2)") from e
        if not token:
            raise SystemExit("no API token: pass --token or set BIOHUB_TOKEN")
        self.model = model
        self.client = esmfold2_client(model=model, url=url, token=token)
        self.config = FoldingConfig(**(config_kwargs or {}))

    def fold(self, sequence, msa=None) -> FoldResult:
        # One ProteinInput chain "A"; msa=None -> single-sequence, else an esm MSA object.
        from esm.utils.structure.input_builder import ProteinInput, StructurePredictionInput
        prot = ProteinInput(id="A", sequence=sequence, msa=msa)
        inp = StructurePredictionInput(sequences=[prot])
        result = self.client.fold_all_atom(inp, config=self.config)
        return _extract_result(sequence, result)


class LocalBackend:
    """Local-GPU inference -- planned as a SEPARATE script, not wired here yet.

    The biohub *tutorial/API* is remote-only (hosted esmfold2_client), but weights ARE
    downloadable from HuggingFace (biohub/ESMFold2-Fast) for local GPU inference:
        from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
        model = ESMFold2Model.from_pretrained("biohub/ESMFold2-Fast").cuda().eval()
        out = model.infer_protein(sequence, num_loops=..., num_sampling_steps=...)
    That path (likely the primary one: a GPU VM with lots of disk) will live in a
    dedicated local script -- see the "Local inference" notes in the issue-#6 notebook.
    For now this backend is disabled; use --backend api.
    """

    def __init__(self, model=DEFAULT_MODEL, device="cuda", config_kwargs=None):
        raise SystemExit(
            "local backend is not wired up here yet -- it's planned as a separate script "
            "(transformers ESMFold2Model.from_pretrained('biohub/ESMFold2-Fast'), see the "
            "issue-#6 notebook). Use --backend api for now.")

    def fold(self, sequence, msa=None):               # pragma: no cover
        raise NotImplementedError


def make_backend(kind, args):
    # FoldingConfig knobs; None entries are omitted so the model's own default applies.
    cfg = {"include_pae": not args.no_pae}
    if args.num_loops is not None:          cfg["num_loops"] = args.num_loops
    if args.num_sampling_steps is not None: cfg["num_sampling_steps"] = args.num_sampling_steps
    if args.lm_mask_pct is not None:        cfg["lm_mask_pct"] = args.lm_mask_pct
    if args.msa_max_depth is not None:      cfg["msa_max_depth"] = args.msa_max_depth
    if kind == "api":
        token = args.token or os.environ.get("BIOHUB_TOKEN")
        return ApiBackend(model=args.model, url=args.url, token=token, config_kwargs=cfg)
    if kind == "local":
        return LocalBackend(model=args.model, device=args.device, config_kwargs=cfg)
    raise SystemExit(f"unknown backend {kind!r}")


# --- MSA arm (petadex MSA+ESM) --- TODO, intentionally unimplemented -------
def build_msa(sequence, centroid_id=None):
    """Return an esm.utils.msa.MSA for the `msa` arm (query first), or None for single-seq.

    Mechanism (now confirmed from the tutorial, Examples 3-4):
        from esm.utils.msa import MSA
        msa = MSA.from_a3m(path=<a3m>, remove_insertions=True, max_sequences=1024)
        #   or MSA.from_sequences([query, homolog1, ...])   # query MUST be first
        # then it flows straight through: ProteinInput(id="A", sequence=seq, msa=msa)
        # (complexes pair rows by a `key=<taxid>` header token; N/A for single centroids)

    TODO: the petadex MSA *source* (precomputed a3m per centroid vs. an on-the-fly
    mmseqs/jackhmmer search of the petadex DB) is still TBD, so --arm msa stays disabled
    until that is wired up -- return an MSA object here and delete this raise when ready.
    """
    raise NotImplementedError(
        "the 'msa' (petadex MSA+ESM) arm is not implemented yet — the MSA source is "
        "still TBD. Run with --arm single for now; return an esm.utils.msa.MSA from "
        "build_msa() once the petadex MSAs are available.")


# --- CSV input -------------------------------------------------------------
def _pick_col(fieldnames, candidates):
    lower = {c.lower().strip(): c for c in fieldnames}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    return None


def read_centroids(csv_path, seq_col=None, id_col=None):
    """Yield (centroid_id, clean_sequence) from the centroid CSV.

    Sequence is upper-cased, stripped of whitespace and a trailing stop '*'. Columns
    are auto-detected from the header unless --seq-col / --id-col pin them. Duplicate
    ids get a numeric suffix so output filenames never collide.
    """
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise SystemExit(f"{csv_path}: empty or headerless CSV")
        sc = seq_col or _pick_col(reader.fieldnames, SEQ_COLS)
        if sc is None:
            if len(reader.fieldnames) == 1:           # single-column file: it's the sequence
                sc = reader.fieldnames[0]
            else:
                raise SystemExit(
                    f"could not find a sequence column in {reader.fieldnames}; "
                    f"pass --seq-col")
        ic = id_col or _pick_col(reader.fieldnames, ID_COLS)
        seen = {}
        for n, row in enumerate(reader):
            raw = (row.get(sc) or "").strip().upper().rstrip("*")
            seq = "".join(raw.split())
            if not seq:
                continue
            bad = set(seq) - VALID_AA
            if bad:
                tqdm.write(f"  WARN row {n}: skipping, non-AA chars {sorted(bad)}")
                continue
            cid = (row.get(ic).strip() if ic and row.get(ic) else f"row{n}")
            cid = _safe_id(cid)
            if cid in seen:                           # de-collide duplicate ids
                seen[cid] += 1
                cid = f"{cid}_{seen[cid]}"
            else:
                seen[cid] = 0
            yield cid, seq


def _safe_id(s):
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:200]


def prot_hash(seq):
    return hashlib.md5(seq.encode()).hexdigest()      # matches the Atlas content-address


# --- output sinks ----------------------------------------------------------
class Sink:
    """Writes structures/<id>.cif and arrays/<id>.npz locally and/or to S3."""

    def __init__(self, outdir, s3_dest=None):
        self.outdir = outdir
        self.struct_dir = os.path.join(outdir, "structures")
        self.array_dir = os.path.join(outdir, "arrays")
        os.makedirs(self.struct_dir, exist_ok=True)
        os.makedirs(self.array_dir, exist_ok=True)
        self.s3 = self.bucket = self.prefix = None
        if s3_dest:
            import boto3
            self.s3 = boto3.client("s3")              # signed; EC2 IAM role writes to petadex
            body = s3_dest[5:]
            self.bucket, _, self.prefix = body.partition("/")
            self.prefix = self.prefix.strip("/")

    def _put_s3(self, subdir, name, data, ctype):
        key = "/".join(p for p in (self.prefix, subdir, name) if p)
        body = data.encode() if isinstance(data, str) else data
        self.s3.put_object(Bucket=self.bucket, Key=key, Body=body, ContentType=ctype)

    def seed_metrics(self, local_path, name="metrics.csv"):
        # Resuming on a fresh box with no local ledger: pull the S3 copy down first so
        # already-folded centroids are skipped. Returns True if a copy was fetched.
        if not self.s3 or os.path.exists(local_path):
            return False
        key = "/".join(p for p in (self.prefix, name) if p)
        try:
            obj = self.s3.get_object(Bucket=self.bucket, Key=key)
        except Exception:
            return False
        with open(local_path, "wb") as f:
            f.write(obj["Body"].read())
        return True

    def mirror_metrics(self, local_path, name="metrics.csv"):
        # Mirror the whole ledger to <prefix>/metrics.csv so it (and resume state) survive an
        # ephemeral VM. Re-uploads the whole file, so call it on an interval + at exit, not per row.
        if not self.s3:
            return
        key = "/".join(p for p in (self.prefix, name) if p)
        with open(local_path, "rb") as f:
            self.s3.put_object(Bucket=self.bucket, Key=key, Body=f.read(), ContentType="text/csv")

    def write(self, cid, fr: FoldResult):
        sname = f"{cid}.cif"
        with open(os.path.join(self.struct_dir, sname), "w") as f:
            f.write(fr.cif_text)

        import io
        L = int(fr.per_residue_plddt.shape[0])
        arrays = dict(per_residue_plddt=fr.per_residue_plddt.astype(np.float32),
                      residue_index=np.arange(1, L + 1, dtype=np.int32))
        if fr.pae is not None:
            arrays["pae"] = fr.pae.astype(np.float32)
        buf = io.BytesIO(); np.savez_compressed(buf, **arrays)
        blob = buf.getvalue()
        aname = f"{cid}.npz"
        with open(os.path.join(self.array_dir, aname), "wb") as f:
            f.write(blob)

        if self.s3:
            self._put_s3("structures", sname, fr.cif_text, "chemical/x-mmcif")
            self._put_s3("arrays", aname, blob, "application/octet-stream")


# --- metrics ledger (also the resume state) --------------------------------
METRIC_COLS = ["centroid_id", "arm", "backend", "model", "seq_len",
               "protein_hash", "wall_s", "mean_plddt", "ptm", "has_pae",
               "status", "error"]


def load_done(metrics_path):
    """Set of (centroid_id, arm) already completed OK, for resume."""
    done = set()
    if os.path.exists(metrics_path):
        with open(metrics_path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("status") == "ok":
                    done.add((row["centroid_id"], row["arm"]))
    return done


# --- bounded concurrent driver ---------------------------------------------
def _bounded(ex, fn, items, max_inflight):
    """Run fn over items with <= max_inflight futures alive; yield results as they finish.
    (Same flat-memory pattern as the Atlas download script.)"""
    it = iter(items)
    inflight = {ex.submit(fn, x) for x in islice(it, max_inflight)}
    while inflight:
        done, inflight = wait(inflight, return_when=FIRST_COMPLETED)
        for f in done:
            yield f.result()
            nxt = next(it, None)
            if nxt is not None:
                inflight.add(ex.submit(fn, nxt))


def fold_one(backend, sink, cid, seq, arm, retries):
    """Fold a single centroid; return a metrics row dict. Never raises (errors captured)."""
    row = {"centroid_id": cid, "arm": arm, "backend": backend.__class__.__name__,
           "model": getattr(backend, "model", ""), "seq_len": len(seq),
           "protein_hash": prot_hash(seq), "wall_s": "", "mean_plddt": "",
           "ptm": "", "has_pae": "", "status": "error", "error": ""}
    msa = build_msa(seq, cid) if arm == "msa" else None   # raises for the unimplemented arm
    t0 = time.perf_counter()
    last = None
    for attempt in range(retries + 1):
        try:
            fr = backend.fold(seq, msa=msa)
            sink.write(cid, fr)
            row.update(wall_s=round(time.perf_counter() - t0, 3),
                       mean_plddt=round(fr.mean_plddt, 3), ptm=round(fr.ptm, 4),
                       has_pae=int(fr.pae is not None), status="ok")
            return row
        except NotImplementedError:
            raise                                     # msa arm: surface immediately
        except Exception as e:                        # transient API/decode failure -> retry
            last = e
            if attempt < retries:
                time.sleep(2 ** attempt)              # 1s, 2s, 4s backoff
    row.update(wall_s=round(time.perf_counter() - t0, 3), error=f"{type(last).__name__}: {last}")
    return row


def run(args):
    records = list(read_centroids(args.csv, seq_col=args.seq_col, id_col=args.id_col))
    if args.limit:
        records = records[:args.limit]
    print(f"{len(records)} centroids from {args.csv} | arm={args.arm} | backend={args.backend}")

    os.makedirs(args.outdir, exist_ok=True)
    metrics_path = os.path.join(args.outdir, "metrics.csv")
    sink = Sink(args.outdir, s3_dest=args.s3)
    if sink.seed_metrics(metrics_path):               # cross-instance resume from S3
        print(f"resume ledger seeded from s3://{sink.bucket}/{sink.prefix}/metrics.csv")
    done = load_done(metrics_path)
    todo = [(c, s) for c, s in records if (c, args.arm) not in done]
    print(f"{len(done)} already done | {len(todo)} to fold")
    if not todo:
        print("nothing to do."); return

    backend = make_backend(args.backend, args)        # constructs client / loads weights

    # Local GPU inference is serial; the remote API is I/O-bound and parallelises well.
    workers = 1 if args.backend == "local" else max(1, args.workers)

    write_header = not os.path.exists(metrics_path)
    ok = fail = 0
    last_sync = time.time()
    with open(metrics_path, "a", newline="") as mf, ThreadPoolExecutor(max_workers=workers) as ex:
        writer = csv.DictWriter(mf, fieldnames=METRIC_COLS)
        if write_header:
            writer.writeheader()
        task = lambda cs: fold_one(backend, sink, cs[0], cs[1], args.arm, args.retries)
        with tqdm(total=len(todo), desc=f"fold[{args.arm}]", unit="seq") as pbar:
            for row in _bounded(ex, task, todo, workers * 2):
                writer.writerow(row); mf.flush()      # crash-safe ledger, one row per fold
                if row["status"] == "ok":
                    ok += 1
                else:
                    fail += 1
                    tqdm.write(f"  FAIL {row['centroid_id']}: {row['error']}")
                if sink.s3 and time.time() - last_sync > 300:   # mirror ledger to S3 ~every 5 min
                    sink.mirror_metrics(metrics_path); last_sync = time.time()
                pbar.update(1)
                pbar.set_postfix(ok=ok, fail=fail)
    sink.mirror_metrics(metrics_path)                 # final mirror so the ledger is durable off-box

    dest = args.s3 if args.s3 else args.outdir
    print(f"\nDone: {ok} folded, {fail} failed. Structures/arrays -> {dest}")
    print(f"Metrics + resume ledger: {metrics_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Fold 60pid centroids with ESMFold2 and benchmark the run.")
    ap.add_argument("csv", help="input CSV of centroid AA sequences (rows = centroids)")
    ap.add_argument("--backend", choices=["api", "local"], default="api",
                    help="api = remote biohub-hosted ESMFold2 (no local GPU); "
                         "local = disabled (the documented SDK is remote-only)")
    ap.add_argument("--arm", choices=["single", "msa"], default="single",
                    help="single = single-sequence fold; msa = petadex MSA+ESM (TODO, disabled)")
    ap.add_argument("--model", default=DEFAULT_MODEL, help=f"model id (default {DEFAULT_MODEL})")
    ap.add_argument("--url", default=DEFAULT_URL, help="biohub API base URL")
    ap.add_argument("--token", default=None, help="API token (else env BIOHUB_TOKEN)")
    ap.add_argument("--device", default="cuda", help="torch device (local backend only)")
    # FoldingConfig knobs (accuracy vs runtime); None leaves the model's own default.
    ap.add_argument("--num-loops", type=int, default=None,
                    help="structure-refinement loops (SDK default 20; fewer = faster/cheaper)")
    ap.add_argument("--num-sampling-steps", type=int, default=None,
                    help="diffusion sampling steps (SDK default 100)")
    ap.add_argument("--lm-mask-pct", type=float, default=None,
                    help="fraction of residues masked before the PLM (fast default 0.1)")
    ap.add_argument("--msa-max-depth", type=int, default=None,
                    help="MSA rows subsampled per loop (SDK default 1024; msa arm only)")
    ap.add_argument("--no-pae", action="store_true",
                    help="skip PAE (faster; arrays/*.npz then omits pae and has_pae=0)")
    ap.add_argument("--outdir", default="out", help="local output dir (default out/)")
    ap.add_argument("--s3", dest="s3", default=S3_DEST, metavar="s3://bucket/prefix/",
                    help=f"also upload structures/ and arrays/ here (default {S3_DEST}); "
                         f"pass --s3 '' to keep output local-only")
    ap.add_argument("--seq-col", default=None, help="sequence column name (else auto-detect)")
    ap.add_argument("--id-col", default=None, help="centroid-id column name (else auto-detect)")
    ap.add_argument("--workers", type=int, default=8,
                    help="concurrent fold requests for the api backend (default 8)")
    ap.add_argument("--retries", type=int, default=3,
                    help="retries per sequence on transient failure (default 3)")
    ap.add_argument("--limit", type=int, default=None, help="only fold the first N centroids")
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
