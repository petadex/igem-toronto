"""
Local ESMFold2 centroid folding -- ONE ESMC backbone per run (this is the script you run).

Architecture (verified against the HF repos, not the prompt's notes):
    biohub/ESMFold2        -> the TRUNK ONLY (~235 MB). Its config.json carries "esmc_id".
    biohub/ESMC-6B         -> the ESMC language-model backbone, a SEPARATE 6.3B-param repo.
    sequence -> ESMC (2560-dim embeddings) -> ESMFold2 trunk -> all-atom structure

So you do NOT swap weights to change the backbone: you tell the trunk which ESMC repo to pair
with. `ESMFold2Model.from_pretrained("biohub/ESMFold2")` reads esmc_id from config and loads that
ESMC; overriding esmc_id points it at a different one (e.g. a finetune). That is what --esmc does.

    baseline run (now):     python esmfold2_local_predictor.py centroids.fasta --outdir out_base
    informed run (later):   python esmfold2_local_predictor.py centroids.fasta --outdir out_ft \
                                --esmc petadex/ESMC-6B-catalytic --label informed

Each run writes structures + a results.csv for ONE backbone. Compare two runs afterwards with
compare_runs.py. If the esmc_id override turns out not to be honoured by the installed esm build,
fall back to --swap-fallback, which loads the base trunk and overwrites the ESMC weights in place
using esmc_backbone_swap.py (same effect, uglier path).

CAVEATS to carry into every interpretation:
  * pLDDT / pTM / ipTM are the model's SELF-CONFIDENCE, not accuracy. A domain-adapted backbone can
    look more confident on catalytic sequences without being more correct. Real quality needs
    experimental references (see --reference-dir) or, for run-vs-run structural change, compare_runs.py.
  * OFF-DISTRIBUTION RISK: a heavy LoRA domain shift can push ESMC embeddings off the distribution
    the frozen trunk was trained to read, degrading folds even though the backbone is "better" at
    catalytic MLM. If the finetuned run's pLDDT drops across the board, that is the likely cause.

Install:  pip install torch transformers "esm @ git+https://github.com/Biohub/esm.git@main" numpy boto3
Example (no GPU/weights, exercises the whole pipeline with a fake model):
          python esmfold2_local_predictor.py centroids.fasta --mock-fold --outdir out_mock

AWS / S3:  structures/, arrays/ and results.csv upload to --s3 (default the petadex sink); the
    results.csv ledger is mirrored to S3 on an interval + at exit and re-seeded from S3 on a fresh
    box, so an ephemeral/spot VM RESUMES cleanly (already-folded ids are skipped). --s3 '' keeps
    output local-only. Confirm weights + measure the dominant (download+load) cost first with:
          python esmfold2_local_predictor.py centroids.fasta --load-only     # loads model, folds nothing
    Validate the S3 path with NO GPU by combining --mock-fold with a real --s3 bucket.
"""
import os
import sys
import time
import argparse
import statistics
from dataclasses import dataclass
from typing import Optional

import numpy as np

ESMFOLD2_REPOS = {"full": "biohub/ESMFold2", "fast": "biohub/ESMFold2-Fast"}
ESMC_6B_REPO = "biohub/ESMC-6B"          # the trunk's default backbone (config esmc_id)
ESMC_6B_HIDDEN = 2560                    # lm_d_model in the trunk config; the dimension it reads
REFERENCE_SPEED = {"full": 15.8, "fast": 9.4}   # H100, 1024 aa -- sanity check for wall times
S3_DEST = "s3://petadex-protein-structures/esmfold2-centroids/"   # canonical sink (issue #6)
MIRROR_EVERY_S = 300                     # re-upload results.csv to S3 at most this often (spot-safe)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
ESMC_CONTEXT = 2048                      # tokens; longer chains get truncated by the backbone


# ===========================================================================
# Sequence input -- FASTA (>headers) OR plain one-seq-per-line (.seqs.txt)
# ===========================================================================
# The petadex ORF/centroid export is one RAW amino-acid sequence per line, no headers, sometimes a
# trailing '*' (stop). That is NOT FASTA, so we auto-detect: a file whose first non-blank line
# starts with '>' is parsed as FASTA; otherwise every non-blank line is its own sequence, with ids
# synthesised from the file basename (e.g. cluster1.seqs.txt -> cluster1_1, cluster1_2, ...).
def _clean_seq(raw):
    return "".join(raw.upper().rstrip("*").split())


def _warn_bad_chars(cid, seq):
    bad = set(seq) - VALID_AA
    if bad:
        print(f"  WARN {cid}: non-AA chars {sorted(bad)} present (kept verbatim)")


def _sniff_format(path):
    with open(path) as fh:
        for line in fh:
            if line.strip():
                return "fasta" if line.lstrip().startswith(">") else "lines"
    return "lines"


def _basename_stem(path):
    stem = os.path.splitext(os.path.basename(path))[0]
    if stem.endswith(".seqs"):        # cluster1.seqs.txt -> cluster1
        stem = stem[:-5]
    return _safe_id(stem) or "seq"


def read_sequences(path, fmt="auto"):
    """Yield (id, clean_sequence). Auto-detects FASTA vs one-seq-per-line; override with fmt."""
    if fmt == "auto":
        fmt = _sniff_format(path)
    if fmt == "fasta":
        yield from _read_fasta(path)
    elif fmt == "lines":
        yield from _read_lines(path)
    else:
        raise SystemExit(f"unknown input format {fmt!r} (use auto|fasta|lines)")


def _read_lines(path):
    """One raw sequence per line; ids = <basename>_<n> (1-based over non-blank lines)."""
    stem = _basename_stem(path)
    n = 0
    with open(path) as fh:
        for line in fh:
            seq = _clean_seq(line.strip())
            if not seq:
                continue
            n += 1
            cid = _safe_id(f"{stem}_{n}")
            _warn_bad_chars(cid, seq)
            yield cid, seq


def _read_fasta(path):
    seen = {}

    def _emit(header, chunks):
        seq = _clean_seq("".join(chunks))
        if not seq:
            return None
        cid = _safe_id(header.split()[0]) if header.split() else "seq"
        if cid in seen:
            seen[cid] += 1
            cid = f"{cid}_{seen[cid]}"
        else:
            seen[cid] = 0
        _warn_bad_chars(cid, seq)
        return cid, seq

    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    rec = _emit(header, chunks)
                    if rec:
                        yield rec
                header, chunks = line[1:], []
            elif line.strip():
                chunks.append(line.strip())
    if header is not None:
        rec = _emit(header, chunks)
        if rec:
            yield rec


def _safe_id(s):
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:200]


def prot_hash(seq):
    import hashlib
    return hashlib.md5(seq.encode()).hexdigest()


def warn_length(cid, seq):
    if len(seq) > ESMC_CONTEXT:
        print(f"  *** WARN {cid}: {len(seq)} residues > ESMC context {ESMC_CONTEXT}. The backbone "
              f"WILL TRUNCATE this chain; structure past residue {ESMC_CONTEXT} is unreliable.")
        return True
    return False


# ===========================================================================
# Determinism + device
# ===========================================================================
def set_seeds(seed):
    import random
    import torch
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_device(requested):
    import torch
    if requested == "cuda" and not torch.cuda.is_available():
        raise SystemExit(
            "no CUDA device available. ESMC-6B needs an 80 GB-class GPU (bf16 backbone alone is "
            "~12-13 GB, plus trunk + activations). Use --mock-fold for a CPU dry run, or run on a GPU.")
    return torch.device(requested)


def dtype_from_str(name):
    import torch
    return {"bf16": torch.bfloat16, "fp16": torch.float16, "fp32": torch.float32}[name]


def cuda_sync(device):
    if getattr(device, "type", None) != "cuda":       # CPU/mock path stays torch-free
        return
    import torch
    torch.cuda.synchronize(device)


# ===========================================================================
# Model loading -- config-driven backbone selection (primary) or swap fallback
# ===========================================================================
def load_esmfold2(variant, esmc_id, device, dtype, swap_fallback=False, verbose=True):
    """Load the ESMFold2 trunk paired with the chosen ESMC backbone.

    PRIMARY path: pass esmc_id as a from_pretrained kwarg so the trunk pairs with that ESMC repo.
    We then print the esmc_id actually in effect so you can confirm the override took.

    FALLBACK path (swap_fallback=True): load the trunk with its DEFAULT ESMC, then overwrite the
    ESMC submodule's weights from `esmc_id` via esmc_backbone_swap.py. Same result, used only if the
    esmc_id override is not honoured by the installed build.
    """
    repo = ESMFOLD2_REPOS[variant]
    try:
        from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    except Exception as e:                            # noqa: BLE001
        raise SystemExit(
            "could not import ESMFold2Model. Install the ESMFold2-aware build:\n"
            "    pip install transformers \"esm @ git+https://github.com/Biohub/esm.git@main\"\n"
            f"underlying error: {type(e).__name__}: {e}") from e

    if swap_fallback and esmc_id:
        if verbose:
            print(f"loading {repo} with DEFAULT backbone, then swapping in {esmc_id} (fallback) ...")
        model = ESMFold2Model.from_pretrained(repo, torch_dtype=dtype).to(device).eval()
        try:
            from esmc_backbone_swap import find_esmc_backbone, swap_in_backbone
        except Exception as e:                        # noqa: BLE001
            raise SystemExit(f"--swap-fallback needs esmc_backbone_swap.py importable: {e}") from e
        _, backbone = find_esmc_backbone(model, verbose=verbose)
        swap_in_backbone(backbone, esmc_id, verbose=verbose)
        return model

    kwargs = {"torch_dtype": dtype}
    if esmc_id:
        kwargs["esmc_id"] = esmc_id                   # config override: which ESMC repo to pair with
    if verbose:
        print(f"loading {repo} (variant={variant}, dtype={dtype}, "
              f"esmc={esmc_id or 'config default'}) ...")
    try:
        model = ESMFold2Model.from_pretrained(repo, **kwargs)
    except Exception as e:                            # noqa: BLE001
        raise SystemExit(f"failed to load {repo}: {type(e).__name__}: {e}\n"
                         "Verify the repo id, ESMC access, and free disk.")
    model = model.to(device).eval()

    got = getattr(getattr(model, "config", None), "esmc_id", None)
    if verbose:
        print(f"  esmc_id in effect: {got}")
    if esmc_id and got not in (esmc_id, None) and got != esmc_id:
        print(f"  *** WARN: requested esmc={esmc_id} but config.esmc_id={got}. The override may not "
              f"have been honoured -- confirm the loaded backbone, or rerun with --swap-fallback.")
    return model


# ===========================================================================
# Folding
# ===========================================================================
def _make_fold_fn(model, seq, params, seed):
    """Zero-arg closure that folds `seq` once. Documented builder path first, then fallbacks."""
    if getattr(model, "_is_mock", False):
        return lambda: mock_fold_result(model, seq, seed)

    try:
        from esm.models.esmfold2 import ProteinInput, StructurePredictionInput, ESMFold2InputBuilder
        builder = ESMFold2InputBuilder()
        spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=seq)])
        return lambda: builder.fold(model, spi, seed=seed, **params)
    except Exception:                                 # noqa: BLE001
        pass
    try:
        from esm.utils.structure.input_builder import ProteinInput, StructurePredictionInput
        from esm.models.esmfold2 import ESMFold2InputBuilder
        builder = ESMFold2InputBuilder()
        spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=seq)])
        return lambda: builder.fold(model, spi, seed=seed, **params)
    except Exception:                                 # noqa: BLE001
        pass
    if hasattr(model, "infer_protein"):
        return lambda: model.infer_protein(seq, seed=seed, **params)
    raise SystemExit(
        "could not construct an ESMFold2 folding call. Expected entry points were not importable "
        "(esm.models.esmfold2.ESMFold2InputBuilder / model.infer_protein). Read the installed esm "
        "source and wire _make_fold_fn accordingly.")


def _np(x):
    if x is None:
        return None
    if hasattr(x, "detach"):
        x = x.detach().cpu().float().numpy()
    return np.asarray(x)


@dataclass
class FoldOut:
    cif_text: str
    plddt: np.ndarray
    ptm: float
    iptm: Optional[float]
    pae: Optional[np.ndarray]

    @property
    def mean_plddt(self):
        return float(np.mean(self.plddt)) if self.plddt.size else float("nan")


def extract_fold(result):
    if isinstance(result, Exception):
        raise result
    cif = result.complex.to_mmcif()
    plddt = _np(result.plddt).astype(np.float32).reshape(-1)
    if plddt.size and float(np.nanmax(plddt)) <= 1.0:
        plddt = plddt * 100.0
    ptm = float(_np(result.ptm))
    iptm = getattr(result, "iptm", None)
    iptm = float(_np(iptm)) if iptm is not None else None
    pae = getattr(result, "pae", None)
    if pae is not None:
        pae = _np(pae).astype(np.float32)
        if pae.ndim == 3 and pae.shape[0] == 1:
            pae = pae[0]
    return FoldOut(cif, plddt, ptm, iptm, pae)


def timed_fold(fold_fn, device, repeats):
    """Run `fold_fn` `repeats` times, syncing CUDA around each; return (result, median_seconds)."""
    times, result = [], None
    for _ in range(repeats):
        cuda_sync(device)
        t0 = time.perf_counter()
        result = fold_fn()
        cuda_sync(device)
        times.append(time.perf_counter() - t0)
    return result, statistics.median(times)


# ===========================================================================
# Mock fold (--mock-fold): whole pipeline on CPU, no weights/GPU/downloads
# ===========================================================================
_THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
          "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
          "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def build_mock_model(device, label="base", hidden=128):
    """Tiny fake model whose weights are seeded from `label`, so different --esmc labels give
    different (fake) structures + confidence -- enough to exercise the pipeline and compare_runs.py."""
    import hashlib
    import torch
    import torch.nn as nn
    g = torch.Generator().manual_seed(int(hashlib.md5(label.encode()).hexdigest(), 16) % (2 ** 31))

    class Mock(nn.Module):
        def __init__(self):
            super().__init__()
            self.language_model = nn.Sequential(*[nn.Linear(hidden, hidden) for _ in range(3)])
            self.trunk = nn.Linear(hidden, hidden // 2)
    m = Mock()
    with torch.no_grad():
        for p in m.parameters():
            p.copy_(torch.randn(p.shape, generator=g))
    m = m.to(device).eval()
    m._is_mock = True
    m._label = label
    return m


def _minimal_cif(seq, coords, plddt):
    lines = ["data_mock", "#", "loop_",
             "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
             "_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
             "_atom_site.label_seq_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
             "_atom_site.Cartn_z", "_atom_site.B_iso_or_equiv"]
    for i, aa in enumerate(seq, start=1):
        x, y, z = coords[i - 1]
        lines.append(f"ATOM {i} C CA {_THREE.get(aa, 'GLY')} A {i} "
                     f"{x:.3f} {y:.3f} {z:.3f} {float(plddt[i-1]):.2f}")
    lines.append("#")
    return "\n".join(lines) + "\n"


def mock_fold_result(model, seq, seed):
    """Fake fold: confidence AND coordinates depend on the (label-seeded) backbone weights, so
    different backbones produce measurably different structures."""
    import hashlib
    import torch
    with torch.no_grad():
        wsum = float(sum(p.float().abs().sum().item() for p in model.language_model.parameters()))
    rng = np.random.default_rng(int(hashlib.md5(f"{seq}|{round(wsum, 2)}|{seed}".encode()).hexdigest(), 16) % (2 ** 32))
    L = len(seq)
    plddt = np.clip(60 + 30 * rng.random(L), 0, 100).astype(np.float32)
    # a wiggly CA trace: ~3.8 A steps perturbed by the backbone-seeded rng
    steps = np.array([3.8, 0.0, 0.0]) + 0.6 * rng.standard_normal((L, 3))
    coords = np.cumsum(steps, axis=0)

    class _Cx:
        def to_mmcif(self_):
            return _minimal_cif(seq, coords, plddt)

    class _R:
        pass
    r = _R()
    r.complex, r.plddt, r.ptm, r.iptm, r.pae = _Cx(), plddt, float(np.clip(0.5 + 0.4 * rng.random(), 0, 1)), None, None
    return r


# ===========================================================================
# Output sink + results ledger
# ===========================================================================
class Sink:
    """Writes structures/<id>.cif and arrays/<id>.npz locally and (optionally) to S3.

    Also mirrors/seeds the results.csv ledger so a benchmark survives an ephemeral
    (spot) VM: seed_results() pulls any existing ledger down on a fresh box so finished
    folds are skipped, mirror_results() re-uploads the whole ledger on an interval + at
    exit. S3 auth is whatever boto3 finds -- on EC2 that's the attached IAM role.
    """

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
            body = s3_dest[5:]                         # strip "s3://"
            self.bucket, _, self.prefix = body.partition("/")
            self.prefix = self.prefix.strip("/")

    def _put_s3(self, subdir, name, data, ctype):
        key = "/".join(p for p in (self.prefix, subdir, name) if p)
        body = data.encode() if isinstance(data, str) else data
        self.s3.put_object(Bucket=self.bucket, Key=key, Body=body, ContentType=ctype)

    def seed_results(self, local_path, name=None):
        """Fresh box, no local ledger: pull the S3 copy so done folds are skipped. The S3 leaf
        defaults to the local basename, so a sharded run seeds only its OWN ledger (no cross-shard
        clobber -- each shard owns a disjoint set of ORFids)."""
        if not self.s3 or os.path.exists(local_path):
            return False
        name = name or os.path.basename(local_path)
        key = "/".join(p for p in (self.prefix, name) if p)
        try:
            obj = self.s3.get_object(Bucket=self.bucket, Key=key)
        except Exception:                             # noqa: BLE001  (no prior run -> nothing to seed)
            return False
        with open(local_path, "wb") as f:
            f.write(obj["Body"].read())
        return True

    def mirror_results(self, local_path, name=None):
        """Mirror the whole ledger to <prefix>/<basename> (interval + at exit, not per row).
        Per-shard basename keeps fleet workers from overwriting each other's ledger."""
        if not self.s3:
            return
        name = name or os.path.basename(local_path)
        key = "/".join(p for p in (self.prefix, name) if p)
        with open(local_path, "rb") as f:
            self.s3.put_object(Bucket=self.bucket, Key=key, Body=f.read(), ContentType="text/csv")

    def write(self, cid, fo: FoldOut):
        cif_name = f"{cid}.cif"
        with open(os.path.join(self.struct_dir, cif_name), "w") as f:
            f.write(fo.cif_text)
        L = int(fo.plddt.shape[0])
        arrays = dict(per_residue_plddt=fo.plddt.astype(np.float32),
                      residue_index=np.arange(1, L + 1, dtype=np.int32))
        if fo.pae is not None:
            arrays["pae"] = fo.pae.astype(np.float32)
        import io
        buf = io.BytesIO(); np.savez_compressed(buf, **arrays); blob = buf.getvalue()
        npz_name = f"{cid}.npz"
        with open(os.path.join(self.array_dir, npz_name), "wb") as f:
            f.write(blob)
        if self.s3:
            self._put_s3("structures", cif_name, fo.cif_text, "chemical/x-mmcif")
            self._put_s3("arrays", npz_name, blob, "application/octet-stream")


def load_done(results_path):
    """Set of centroid_ids already folded OK (single-backbone run -> keyed by id), for resume."""
    import csv
    done = set()
    if os.path.exists(results_path):
        with open(results_path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("status") == "ok":
                    done.add(row["centroid_id"])
    return done


RESULT_COLS = ["centroid_id", "run_label", "esmfold2", "esmc", "seq_len", "protein_hash",
               "truncated", "wall_s_median", "residues_per_s", "mean_plddt", "ptm", "iptm",
               "status", "error"]


@dataclass
class Row:
    centroid_id: str
    run_label: str
    esmfold2: str
    esmc: str
    seq_len: int
    protein_hash: str
    truncated: int = 0
    wall_s_median: object = ""
    residues_per_s: object = ""
    mean_plddt: object = ""
    ptm: object = ""
    iptm: object = ""
    status: str = "error"
    error: str = ""

    def as_dict(self):
        return {k: getattr(self, k) for k in RESULT_COLS}


def write_results(path, rows):
    import csv
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=RESULT_COLS)
        w.writeheader()
        for r in rows:
            w.writerow(r.as_dict())


def _coerce_num(v):
    """CSV cells are strings; turn numeric ones back into float so the summary can median them."""
    if v is None or v == "":
        return ""
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def read_result_rows(path):
    """Reconstruct Row objects from an existing results.csv (for the end-of-run summary on resume)."""
    import csv
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, newline="") as fh:
        for d in csv.DictReader(fh):
            rows.append(Row(
                centroid_id=d.get("centroid_id", ""), run_label=d.get("run_label", ""),
                esmfold2=d.get("esmfold2", ""), esmc=d.get("esmc", ""),
                seq_len=int(_coerce_num(d.get("seq_len")) or 0),
                protein_hash=d.get("protein_hash", ""),
                truncated=int(_coerce_num(d.get("truncated")) or 0),
                wall_s_median=_coerce_num(d.get("wall_s_median")),
                residues_per_s=_coerce_num(d.get("residues_per_s")),
                mean_plddt=_coerce_num(d.get("mean_plddt")), ptm=_coerce_num(d.get("ptm")),
                iptm=_coerce_num(d.get("iptm")), status=d.get("status", "error"),
                error=d.get("error", "")))
    return rows


def _median_of(rows, attr):
    vals = [getattr(r, attr) for r in rows
            if r.status == "ok" and isinstance(getattr(r, attr), (int, float))]
    return statistics.median(vals) if vals else float("nan")


def print_summary(rows, variant, label):
    ok = sum(1 for r in rows if r.status == "ok")
    print("\n" + "=" * 70)
    print(f"SUMMARY  run='{label}'  ({ok}/{len(rows)} folded ok)  medians:")
    print("=" * 70)
    for lbl, attr, fmt in [("wall_s (median)", "wall_s_median", "{:.3f}"),
                           ("residues/sec", "residues_per_s", "{:.1f}"),
                           ("mean pLDDT", "mean_plddt", "{:.2f}"),
                           ("pTM", "ptm", "{:.3f}")]:
        v = _median_of(rows, attr)
        print(f"  {lbl:<20} {fmt.format(v) if v == v else 'nan'}")
    print(f"\nReference speed (H100, 1024 aa): ESMFold2-{variant} ~= {REFERENCE_SPEED[variant]:.1f}s.")
    print("NOTE: pLDDT/pTM are SELF-CONFIDENCE, not accuracy. Compare two runs (baseline vs finetune)")
    print("      with compare_runs.py; for real accuracy score against experimental refs.")
    print("=" * 70)


# ===========================================================================
# Reference scoring (vs EXPERIMENTAL structures) -- integration seam
# ===========================================================================
def score_against_reference(pred_cif_path, reference_dir, cid):
    """Score a prediction against an experimental reference (TM-score/lDDT/RMSD). Integration seam:
    wire US-align/TMalign or biotite here. Missing reference/tooling returns None -- never faked.
    NB: for run-vs-run (baseline vs finetune) structural change, use compare_runs.py instead."""
    for ext in (".cif", ".pdb", ".mmcif"):
        cand = os.path.join(reference_dir, cid + ext)
        if os.path.exists(cand):
            print(f"  [ref] found {cand} for {cid} -- scoring not wired (seam in score_against_reference).")
            return None
    return None


# ===========================================================================
# Per-run driver
# ===========================================================================
def fold_run(model, centroids, sink, results_path, done, label, variant, esmc,
             device, params, seed, warmup, repeats):
    """Fold centroids, writing one results.csv row per fold (crash-safe) and mirroring the
    ledger to S3 on an interval + at exit, so a reclaimed spot VM resumes cleanly.
    `done` = centroid_ids already folded OK (skipped). Returns this session's new rows."""
    import csv
    todo = [(c, s) for c, s in centroids if c not in done]
    print(f"\n--- folding {len(todo)} centroids ({len(done)} already done)  run='{label}'  "
          f"esmc={esmc}  (warmup={warmup}, repeats={repeats}, seed={seed}) ---")
    if not todo:
        return []

    if warmup:
        wcid, wseq = todo[0]
        print(f"  warm-up fold on {wcid} ({len(wseq)} aa, not timed) ...")
        try:
            timed_fold(_make_fold_fn(model, wseq, params, seed), device, 1)
        except Exception as e:                        # noqa: BLE001
            print(f"  warm-up failed ({type(e).__name__}: {e}) -- continuing.")

    rows = []
    write_header = not os.path.exists(results_path) or os.path.getsize(results_path) == 0
    last_mirror = time.time()
    with open(results_path, "a", newline="") as mf:
        writer = csv.DictWriter(mf, fieldnames=RESULT_COLS)
        if write_header:
            writer.writeheader()
        for cid, seq in todo:
            row = Row(centroid_id=cid, run_label=label, esmfold2=variant, esmc=esmc, seq_len=len(seq),
                      protein_hash=prot_hash(seq), truncated=int(warn_length(cid, seq)))
            try:
                result, med = timed_fold(_make_fold_fn(model, seq, params, seed), device, repeats)
                fo = extract_fold(result)
                sink.write(cid, fo)
                row.wall_s_median = round(med, 3)
                row.residues_per_s = round(len(seq) / med, 2) if med > 0 else ""
                row.mean_plddt = round(fo.mean_plddt, 3)
                row.ptm = round(fo.ptm, 4)
                row.iptm = round(fo.iptm, 4) if fo.iptm is not None else ""
                row.status = "ok"
                print(f"  {cid:<24} L={len(seq):<5} {med:7.2f}s  pLDDT={fo.mean_plddt:6.2f}  pTM={fo.ptm:.4f}")
            except Exception as e:                    # noqa: BLE001
                row.error = f"{type(e).__name__}: {e}"
                print(f"  FAIL {cid}: {row.error}")
            writer.writerow(row.as_dict()); mf.flush()   # crash-safe: one row per fold
            rows.append(row)
            if sink.s3 and time.time() - last_mirror > MIRROR_EVERY_S:
                sink.mirror_results(results_path); last_mirror = time.time()
    sink.mirror_results(results_path)                 # durable off-box before we exit
    return rows


# ===========================================================================
# Sharding -- one image fans out across a fleet unchanged (CLI --shard or $ESMFOLD2_SHARD)
# ===========================================================================
def parse_shard(spec):
    """'K/N' (1-based) -> (k, n); None for empty input."""
    if not spec:
        return None
    try:
        ks, ns = str(spec).split("/")
        k, n = int(ks), int(ns)
    except Exception:                                 # noqa: BLE001
        raise SystemExit(f"--shard must be 'K/N' (1-based), got {spec!r}")
    if not (1 <= k <= n):
        raise SystemExit(f"--shard 'K/N' needs 1<=K<=N, got {k}/{n}")
    return k, n


def select_shard(centroids, k, n):
    """Sort by sequence length, return the contiguous band K of N. Length-banding (NOT a random
    slice) keeps each worker inside one/two padding buckets, so the compile cache stays warm and
    throughput is uniform. Ordering is deterministic (length, then id) so shards never overlap or
    gap regardless of input order."""
    ordered = sorted(centroids, key=lambda cs: (len(cs[1]), cs[0]))
    bounds = [round(i * len(ordered) / n) for i in range(n + 1)]
    return ordered[bounds[k - 1]:bounds[k]]


def run(args):
    import torch

    if args.hf_cache:                                 # point HF at a persistent/pre-warmed cache
        os.environ["HF_HOME"] = args.hf_cache
        os.environ["HF_HUB_CACHE"] = os.path.join(args.hf_cache, "hub")
        print(f"HF cache -> {args.hf_cache}")

    centroids = list(read_sequences(args.fasta, fmt=args.in_format))
    detected = _sniff_format(args.fasta) if args.in_format == "auto" else args.in_format
    print(f"{len(centroids)} sequences from {args.fasta} (format: {detected})")

    shard = parse_shard(args.shard if args.shard is not None else os.environ.get("ESMFOLD2_SHARD"))
    tag = ""
    if shard:
        k, n = shard
        before = len(centroids)
        centroids = select_shard(centroids, k, n)
        tag = f".shard{k}of{n}"
        lens = [len(s) for _, s in centroids]
        band = f"{min(lens)}-{max(lens)}aa" if lens else "empty"
        print(f"shard {k}/{n}: {len(centroids)}/{before} sequences, length band {band}")
    if args.limit:
        centroids = centroids[:args.limit]
    if not centroids:
        raise SystemExit(f"no usable sequences after shard/limit on {args.fasta}")

    label = args.label or (_safe_id(args.esmc.split('/')[-1]) if args.esmc else "baseline")
    esmc_repr = args.esmc or ESMC_6B_REPO
    device = torch.device("cpu") if args.mock_fold else resolve_device(args.device)
    dtype = dtype_from_str(args.dtype)
    set_seeds(args.seed)

    params = dict(num_loops=args.num_loops, num_sampling_steps=args.num_sampling_steps,
                  num_diffusion_samples=args.num_diffusion_samples)
    print(f"fold params: {params}")
    os.makedirs(args.outdir, exist_ok=True)

    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)    # clean baseline for the peak-VRAM read
    t_load = time.perf_counter()
    if args.mock_fold:
        print("*** MOCK-FOLD: no real weights/GPU. Full pipeline with a tiny label-seeded fake model;\n"
              "    structures/confidence are NOT real -- this validates plumbing only.")
        model = build_mock_model(device, label=esmc_repr)
    else:
        model = load_esmfold2(args.esmfold2, args.esmc, device, dtype,
                              swap_fallback=args.swap_fallback, verbose=True)
    load_s = time.perf_counter() - t_load
    print(f"model load (+download if first run): {load_s:.1f}s")
    if device.type == "cuda":
        vram = torch.cuda.max_memory_allocated(device) / 1e9
        print(f"peak VRAM after load: {vram:.1f} GB")

    if args.load_only:
        print("--load-only: model loaded, no folding. Use this to confirm weights resolve and to\n"
              "             measure download+load cost before committing to a full run.")
        return

    sink = Sink(args.outdir, s3_dest=(args.s3 or None))
    results_name = f"results{tag}.csv"                # per-shard ledger key -> no fleet clobber
    results_path = os.path.join(args.outdir, results_name)
    if sink.seed_results(results_path):               # cross-instance resume from S3 (own shard only)
        print(f"resume ledger seeded from s3://{sink.bucket}/{sink.prefix}/{results_name}")
    done = load_done(results_path)

    rows = fold_run(model, centroids, sink, results_path, done, label, args.esmfold2, esmc_repr,
                    device, params, args.seed, args.warmup, args.repeats)

    if args.reference_dir:
        print(f"\nscoring vs references in {args.reference_dir} ...")
        for r in rows:
            if r.status == "ok":
                score_against_reference(os.path.join(args.outdir, "structures", f"{r.centroid_id}.cif"),
                                        args.reference_dir, r.centroid_id)

    dest = args.s3 if args.s3 else results_path
    print(f"\nresults -> {results_path}   (structures/arrays/ledger -> {dest})")
    print_summary(read_result_rows(results_path), args.esmfold2, label)
    print(f"\nNext: run the finetune later into a different --outdir, then\n"
          f"      python compare_runs.py {args.outdir} <finetune_outdir>")


# ===========================================================================
# Self-test (no GPU/weights)
# ===========================================================================
def self_test():
    import tempfile
    print("SELF-TEST: FASTA parsing, timing, reporting (no GPU/weights).")
    with tempfile.TemporaryDirectory() as td:
        fp = os.path.join(td, "t.fa")
        with open(fp, "w") as f:
            f.write(">c1 d\nMKT AR\nGGP*\n>c1\nMMMM\n>c2\n" + "A" * (ESMC_CONTEXT + 5) + "\n")
        assert _sniff_format(fp) == "fasta"
        recs = list(read_sequences(fp))
        assert [r[0] for r in recs] == ["c1", "c1_1", "c2"], recs
        assert recs[0][1] == "MKTARGGP", recs[0][1]
        assert warn_length("c2", recs[2][1]) is True
    print("  [ok] FASTA parse + de-collision + length warning")

    # plain one-seq-per-line (.seqs.txt), the petadex ORF/centroid export format
    with tempfile.TemporaryDirectory() as td:
        fp = os.path.join(td, "cluster1.seqs.txt")
        with open(fp, "w") as f:
            f.write("MKTAYIAK\nGGSPLLVQ*\n\n  \nqrstvwy*\n")   # blank lines + trailing '*' + lowercase
        assert _sniff_format(fp) == "lines"
        recs = list(read_sequences(fp))
        assert [r[0] for r in recs] == ["cluster1_1", "cluster1_2", "cluster1_3"], recs
        assert recs[1][1] == "GGSPLLVQ" and recs[2][1] == "QRSTVWY", recs
    print("  [ok] one-seq-per-line parse + basename ids + '*'/blank/lowercase handling")

    class _Dev:
        type = "cpu"
    n = {"c": 0}

    def fake():
        n["c"] += 1
        return "r"
    res, med = timed_fold(fake, _Dev(), repeats=3)
    assert res == "r" and n["c"] == 3 and med >= 0
    print("  [ok] timed_fold: N repeats + median")

    rows = [Row("c1", "baseline", "full", ESMC_6B_REPO, 300, "h", 0, 12.0, 25.0, 88.0, 0.85, "", "ok", "")]
    with tempfile.TemporaryDirectory() as td:
        rp = os.path.join(td, "r.csv")
        write_results(rp, rows)
        assert os.path.exists(rp)
    print_summary(rows, "full", "baseline")
    print("  [ok] results.csv + summary")
    print("\nSELF-TEST PASSED")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Local ESMFold2 folding -- one ESMC backbone per run (see compare_runs.py to diff runs).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", nargs="?", metavar="input",
                    help="input sequences: FASTA (>headers) or one raw seq per line (.seqs.txt)")
    ap.add_argument("--in-format", choices=["auto", "fasta", "lines"], default="auto",
                    help="input format (default auto: '>' first line = fasta, else one-seq-per-line)")
    ap.add_argument("--outdir", default="out", help="output dir (default out/)")
    ap.add_argument("--label", default=None,
                    help="run label recorded in results.csv (default: 'baseline' or the esmc basename)")

    ap.add_argument("--esmfold2", choices=["full", "fast"], default="full",
                    help="full=biohub/ESMFold2 (48 trunk layers); fast=ESMFold2-Fast (24)")
    ap.add_argument("--esmc", default=None,
                    help="ESMC backbone repo/path to pair with the trunk (esmc_id override). "
                         "Omit for the trunk's default (biohub/ESMC-6B). Point at the finetune for "
                         "the informed run.")
    ap.add_argument("--swap-fallback", action="store_true",
                    help="if the esmc_id override is not honoured: load the base trunk and overwrite "
                         "the ESMC weights from --esmc via esmc_backbone_swap.py")
    ap.add_argument("--device", default="cuda", help="torch device (default cuda)")
    ap.add_argument("--dtype", choices=["bf16", "fp16", "fp32"], default="bf16", help="compute dtype")

    ap.add_argument("--num-loops", type=int, default=20, help="structure-refinement loops (default 20)")
    ap.add_argument("--num-sampling-steps", type=int, default=100, help="diffusion steps (default 100)")
    ap.add_argument("--num-diffusion-samples", type=int, default=1, help="diffusion samples (default 1)")

    ap.add_argument("--seed", type=int, default=0, help="global seed (fixes diffusion sampling)")
    ap.add_argument("--repeats", type=int, default=1, help="timed repeats per sequence; median reported")
    ap.add_argument("--no-warmup", dest="warmup", action="store_false", help="skip the untimed warm-up")
    ap.add_argument("--reference-dir", default=None,
                    help="dir of experimental structures (<id>.cif/.pdb) for the scoring seam")
    ap.add_argument("--limit", type=int, default=None, help="only fold the first N centroids")
    ap.add_argument("--shard", default=None, metavar="K/N",
                    help="fold only length-band K of N (1-based) so one image fans out across a "
                         "fleet unchanged; sequences are sorted by length and split into N "
                         "contiguous bands (warms one/two padding buckets per worker). Falls back "
                         "to $ESMFOLD2_SHARD (e.g. from instance user-data). Each shard writes its "
                         "own results.shardKofN.csv ledger, so workers never clobber each other.")
    ap.add_argument("--hf-cache", default=None, metavar="DIR",
                    help="set HF_HOME/HF_HUB_CACHE to DIR before loading -- point at a persistent "
                         "EBS volume or a pre-warmed cache so a 6B checkpoint isn't re-downloaded "
                         "on every boot (or bake it into the AMI; see AWS_RUNBOOK.md)")
    ap.add_argument("--s3", dest="s3", default=S3_DEST, metavar="s3://bucket/prefix/",
                    help=f"also upload structures/, arrays/ and results.csv here (default {S3_DEST}); "
                         f"the ledger is mirrored on an interval + at exit so a reclaimed spot VM "
                         f"resumes cleanly. Pass --s3 '' to keep output local-only.")
    ap.add_argument("--load-only", action="store_true",
                    help="load the model, report download+load wall time and peak VRAM, then exit "
                         "(no folding) -- confirm weights resolve + measure the dominant cost first")
    ap.add_argument("--mock-fold", action="store_true",
                    help="run the full pipeline on CPU with a tiny FAKE model (no weights/GPU/download)")
    ap.add_argument("--self-test", action="store_true", help="plumbing checks, then exit")

    args = ap.parse_args()
    if args.self_test:
        sys.exit(self_test())
    if not args.fasta:
        ap.error("a FASTA input is required (or use --self-test)")
    run(args)


if __name__ == "__main__":
    main()
