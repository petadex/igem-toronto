#!/usr/bin/env python
"""Fold protein sequences with ESMFold2 (biohub/ESMFold2 + ESMC-6B backbone).

Self-contained, for the Compute Canada / Vector L40S cluster. One model load, fold a
shard of the input FASTA, write structures + metrics locally. Resume is implicit: any
sequence whose structures/<id>.cif already exists is skipped.

    python fold.py input.fasta --outdir out --shard 3/16 --dtype bf16

Output layout (under --outdir):
    structures/<id>.cif        all-atom mmCIF, per-residue pLDDT in the B-factor column
    metrics/<id>.json          {id, seq_len, params, timing, mean_plddt, ptm, iptm, status}
    progress.shard<K>of<N>.csv one row per attempted fold (append-only, for throughput/ETA)
"""
import argparse
import csv
import gzip
import json
import os
import time


def read_fasta(path):
    """Yield (id, sequence). id = first whitespace token of the header (already >orf<N>)."""
    rid, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(chunks)
                rid = line[1:].split()[0] if line[1:].split() else "seq"
                chunks = []
            else:
                chunks.append(line.strip())
    if rid is not None:
        yield rid, "".join(chunks)


def clean_seq(s):
    return "".join(s.upper().split()).rstrip("*")


# fold cost grows super-linearly with length (~L^1.2 for ESMFold2-Fast) plus a fixed per-fold
# floor, and the set is right-skewed, so equal-COUNT length bands would make the longest-length
# shard ~15x slower than the shortest. Balance shards by predicted cost instead -> equal wall-time.
COST_EXP = 1.2
COST_FLOOR = 100.0    # fixed-overhead term in the same units as L**COST_EXP (~1s at Fast's scale)


def select_shard(records, k, n, cost_exp=COST_EXP, cost_floor=COST_FLOOR):
    """Deterministic cost-balanced partition: sort by (len, id), split so each of the n shards
    carries ~equal total predicted fold cost. Contiguous, disjoint, gapless; keeps each worker
    in a narrow length band (warm compile cache)."""
    import bisect, itertools
    ordered = sorted(records, key=lambda r: (len(r[1]), r[0]))
    cum = list(itertools.accumulate(len(s) ** cost_exp + cost_floor for _, s in ordered))
    total = cum[-1] if cum else 0.0
    lo = bisect.bisect_left(cum, (k - 1) * total / n)
    hi = len(ordered) if k == n else bisect.bisect_left(cum, k * total / n)
    return ordered[lo:hi]


def parse_shard(spec):
    k, n = (int(x) for x in str(spec).split("/"))
    if not (1 <= k <= n):
        raise SystemExit(f"--shard needs 1<=K<=N, got {spec}")
    return k, n


def main():
    ap = argparse.ArgumentParser(description="Fold sequences with ESMFold2.")
    ap.add_argument("input", help="input FASTA (>orf<id> headers)")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--shard", default="1/1", metavar="K/N", help="length-banded shard, 1-based")
    ap.add_argument("--limit", type=int, default=None, help="fold only first N of the shard")
    ap.add_argument("--sample", type=int, default=None,
                    help="fold N sequences evenly spaced across the length range "
                         "(length-stratified, for a representative timing sample)")
    ap.add_argument("--max-len", type=int, default=1280,
                    help="skip sequences longer than this. Measured fp32 OOM ceiling on the "
                         "46 GB L40S is ~1300aa (1285 folds, 1529 OOMs); 1280 keeps 99.75%% of "
                         "the set and avoids the ~25-min allocator stalls a failed fold costs.")
    ap.add_argument("--num-sampling-steps", type=int, default=68)
    ap.add_argument("--num-loops", type=int, default=10)
    # bf16 errors on this build (fp32 tensor meets bf16 weight in the fold path); fp32 is the default.
    ap.add_argument("--dtype", choices=["bf16", "fp16", "fp32"], default="fp32")
    ap.add_argument("--fast", action="store_true",
                    help="use biohub/ESMFold2-Fast (lighter 24-layer trunk) instead of the "
                         "full 48-layer biohub/ESMFold2. Same ESMC-6B backbone, same fold API.")
    # --backend hf is the Biohub transformers fork in ~/petabite (the Aug 2026 calibration).
    # --backend native is esm>=3.4.1's own EsmFold2Model in ~/esmfold2-341: bf16 trunk autocast,
    # chunked L^2 ops and fused kernels. --dtype applies to hf only; native uses --esmc-precision.
    ap.add_argument("--backend", choices=["hf", "native"], default="hf")
    ap.add_argument("--esmc-precision", choices=["bf16", "fp32", "fp8"], default="bf16",
                    help="native only: ESMC-6B backbone precision (fp8 needs H100 + TransformerEngine)")
    ap.add_argument("--kernel", choices=["none", "fused", "cuequivariance"], default="none",
                    help="native only: trunk kernel backend")
    ap.add_argument("--chunk", default="auto",
                    help="native only: token-axis chunk for L^2 ops; 'none' (fastest, most memory), "
                         "an int, or 'auto' = none up to --chunk-above aa, 64 beyond it")
    ap.add_argument("--chunk-above", type=int, default=1000,
                    help="native --chunk auto: sequence length above which chunk=64 is used")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--load-only", action="store_true", help="load model, report load/VRAM, exit")
    # production options (Part II, Step 15). Defaults keep the benchmark behaviour unchanged.
    ap.add_argument("--cost-exp", type=float, default=COST_EXP,
                    help="shard balancing: per-fold cost ~ L**cost_exp + cost_floor. Fit per tier "
                         "(Fast 32/3 >=150aa: ~1.66; Full 32/3 <150aa: ~0.16, i.e. flat)")
    ap.add_argument("--cost-floor", type=float, default=COST_FLOOR)
    ap.add_argument("--gzip", action="store_true", help="write structures/<id>.cif.gz (~4x smaller)")
    ap.add_argument("--subdirs", type=int, default=0,
                    help="spread outputs over N subdirectories (<orf number> %% N), so no single "
                         "directory holds millions of files. 0 = flat layout")
    args = ap.parse_args()

    import torch
    try:  # silence RDKit's per-component CCD-pickle version warning (was 99.9% of the log)
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    from esm.models.esmfold2 import (
        ESMFold2InputBuilder,
        ProteinInput,
        StructurePredictionInput,
    )

    dtype = {"bf16": torch.bfloat16, "fp16": torch.float16, "fp32": torch.float32}[args.dtype]
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise SystemExit("--device cuda but CUDA is not available (need a GPU node).")

    # --- load model once ---
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    t0 = time.perf_counter()
    repo = "biohub/ESMFold2-Fast" if args.fast else "biohub/ESMFold2"
    if args.backend == "hf":
        from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
        model = ESMFold2Model.from_pretrained(repo, torch_dtype=dtype).to(device).eval()
    else:
        from esm.models.esmfold2 import EsmFold2Model
        model = EsmFold2Model.from_pretrained(
            repo, device=device, esmc_precision=args.esmc_precision).eval()
        model.set_kernel_backend(None if args.kernel == "none" else args.kernel)
        if args.chunk != "auto":
            model.set_chunk_size(None if args.chunk == "none" else int(args.chunk))
    print(f"model loaded in {time.perf_counter() - t0:.1f}s ({repo}, backend={args.backend})", flush=True)
    if device.type == "cuda":
        print(f"peak VRAM after load: {torch.cuda.max_memory_allocated(device) / 1e9:.1f} GB", flush=True)
    if args.load_only:
        return

    builder = ESMFold2InputBuilder()

    # --- select this worker's slice ---
    k, n = parse_shard(args.shard)
    records = [(rid, clean_seq(seq)) for rid, seq in read_fasta(args.input) if seq]
    records = select_shard(records, k, n, args.cost_exp, args.cost_floor)
    if args.sample and args.sample < len(records):
        # records are length-sorted, so evenly-spaced indices span the length range
        idx = sorted({round(i * (len(records) - 1) / (args.sample - 1)) for i in range(args.sample)})
        records = [records[j] for j in idx]
    if args.limit:
        records = records[:args.limit]
    if not records:
        raise SystemExit("no sequences after shard/sample/limit")
    tag = f".shard{k}of{n}"
    lens = [len(s) for _, s in records]
    print(f"shard {k}/{n}: {len(records)} sequences, length band {min(lens)}-{max(lens)}aa", flush=True)

    struct_dir = os.path.join(args.outdir, "structures")
    metrics_dir = os.path.join(args.outdir, "metrics")
    os.makedirs(struct_dir, exist_ok=True)
    os.makedirs(metrics_dir, exist_ok=True)
    progress_path = os.path.join(args.outdir, f"progress{tag}.csv")
    progress_cols = ["id", "seq_len", "wall_s", "residues_per_s", "peak_vram_gb", "mean_plddt", "ptm", "iptm", "status"]
    new_progress = not os.path.exists(progress_path)
    pf = open(progress_path, "a", newline="")
    pw = csv.DictWriter(pf, fieldnames=progress_cols)
    if new_progress:
        pw.writeheader()

    params = dict(num_loops=args.num_loops, num_sampling_steps=args.num_sampling_steps)
    print(f"fold params: {params}", flush=True)

    done = skipped = failed = 0
    t_start = time.perf_counter()
    for i, (rid, seq) in enumerate(records, 1):
        sub = ""
        if args.subdirs:
            num = "".join(ch for ch in rid if ch.isdigit())
            sub = f"{int(num) % args.subdirs:03d}" if num else "x"
        os.makedirs(os.path.join(struct_dir, sub), exist_ok=True)
        os.makedirs(os.path.join(metrics_dir, sub), exist_ok=True)
        cif_path = os.path.join(struct_dir, sub, f"{rid}.cif" + (".gz" if args.gzip else ""))
        json_path = os.path.join(metrics_dir, sub, f"{rid}.json")
        if os.path.exists(cif_path) or os.path.exists(cif_path[:-3] if args.gzip else cif_path + ".gz"):
            skipped += 1
            continue
        if len(seq) > args.max_len:
            skipped += 1
            pw.writerow({"id": rid, "seq_len": len(seq), "status": "skipped_too_long"})
            pf.flush()
            print(f"[{i}/{len(records)}] {rid} len={len(seq)} skipped_too_long (>{args.max_len})", flush=True)
            continue
        try:
            if device.type == "cuda":
                torch.cuda.reset_peak_memory_stats(device)
            t = time.perf_counter()
            if args.backend == "native" and args.chunk == "auto":
                model.set_chunk_size(None if len(seq) <= args.chunk_above else 64)
            spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=seq)])
            result = builder.fold(model, spi, seed=args.seed, **params)
            wall = time.perf_counter() - t

            plddt = result.plddt.float().cpu()
            mean_plddt = float(plddt.mean())
            if mean_plddt <= 1.0:            # normalize 0-1 -> 0-100 if needed
                mean_plddt *= 100.0
            ptm = float(result.ptm) if result.ptm is not None else None
            iptm = float(result.iptm) if result.iptm is not None else None

            opener = (lambda p: gzip.open(p, "wt", compresslevel=6)) if args.gzip else (lambda p: open(p, "w"))
            with opener(cif_path + ".tmp") as f:          # write-then-rename: a killed job never
                f.write(result.complex.to_mmcif())        # leaves a truncated file that resume would skip
            metrics = {
                "id": rid, "seq_len": len(seq),
                "num_sampling_steps": args.num_sampling_steps, "num_loops": args.num_loops,
                "dtype": args.dtype, "seed": args.seed, "backend": args.backend,
                "esmc_precision": args.esmc_precision if args.backend == "native" else None,
                "kernel": args.kernel if args.backend == "native" else None,
                "chunk": args.chunk if args.backend == "native" else None,
                "fast": args.fast,
                "peak_vram_gb": round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
                                if device.type == "cuda" else None,
                "wall_s": round(wall, 3), "residues_per_s": round(len(seq) / wall, 2),
                "mean_plddt": round(mean_plddt, 3),
                "ptm": round(ptm, 4) if ptm is not None else None,
                "iptm": round(iptm, 4) if iptm is not None else None,
                "status": "ok",
            }
            with open(json_path, "w") as f:
                json.dump(metrics, f, separators=(",", ":"))
            os.replace(cif_path + ".tmp", cif_path)      # structure appears last: resume keys on it
            pw.writerow({c: metrics.get(c, "") for c in progress_cols})
            pf.flush()
            done += 1
            status = "ok"
        except Exception as e:
            failed += 1
            status = f"error: {type(e).__name__}"
            pw.writerow({"id": rid, "seq_len": len(seq), "status": status})
            pf.flush()
            print(f"  [{rid}] FAILED: {type(e).__name__}: {e}", flush=True)

        elapsed = time.perf_counter() - t_start
        rate = done / elapsed if done else 0
        remaining = len(records) - i
        eta_h = (remaining / rate / 3600) if rate else float("inf")
        print(f"[{i}/{len(records)}] {rid} len={len(seq)} {status} "
              f"| {done} done {skipped} skip {failed} fail "
              f"| {rate * 3600:.0f} folds/h | ETA {eta_h:.1f}h", flush=True)

    pf.close()
    print(f"\nshard {k}/{n} complete: {done} folded, {skipped} skipped, {failed} failed", flush=True)


if __name__ == "__main__":
    main()
