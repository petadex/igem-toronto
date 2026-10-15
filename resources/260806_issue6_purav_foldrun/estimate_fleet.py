#!/usr/bin/env python
"""Extrapolate GPU-hours and storage for the full 90pid fold from a stratified test sample.

Fold time is super-linear in sequence length, and the 90pid centroids span 29..59,275 aa
(median 264). A fixed-throughput estimate from one length is useless. So this fits a power
law wall_s = a * L^b from a LENGTH-STRATIFIED sample (fold.py --sample), then integrates
that over the true length histogram to get total GPU-time.

    # 1. fold a stratified sample spanning the length range
    python fold.py data/90pid_catalytic_orfs.orfid.fa --outdir data/sample --sample 80 --max-len 2048
    # 2. build the length histogram (length<TAB>count), if not already present
    #    seqkit fx2tab -nil INPUT | awk '{c[$NF]++} END{for(l in c) print l"\\t"c[l]}' | sort -n > data/length_hist.tsv
    # 3. estimate
    python estimate_fleet.py --test-dir data/sample --length-hist data/length_hist.tsv
"""
import argparse
import csv
import glob
import math
import os


def dir_bytes(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            total += os.path.getsize(os.path.join(root, f))
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-dir", required=True, help="a completed stratified sample --outdir")
    ap.add_argument("--length-hist", required=True, help="TSV: length<TAB>count over the full input")
    ap.add_argument("--max-len", type=int, default=1280, help="sequences longer than this are skipped (fp32 L40S ceiling)")
    ap.add_argument("--disk-budget-gb", type=float, default=740.0)
    ap.add_argument("--gpu-scenarios", default="8,16,32,64,128")
    args = ap.parse_args()

    import numpy as np

    # --- sample: (seq_len, wall_s) from the stratified fold ---
    rows = []
    for p in glob.glob(os.path.join(args.test_dir, "progress*.csv")):
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r.get("status") == "ok" and r.get("wall_s"):
                    rows.append((int(r["seq_len"]), float(r["wall_s"])))
    if len(rows) < 3:
        raise SystemExit(f"need >=3 completed folds in {args.test_dir}; got {len(rows)}")
    L = np.array([r[0] for r in rows], float)
    W = np.array([r[1] for r in rows], float)
    if L.max() / L.min() < 3:
        print(f"WARNING: sample length range is narrow ({int(L.min())}-{int(L.max())}aa). "
              "Re-fold with fold.py --sample to span the full range, or the fit is unreliable.\n")

    # power-law fit  wall = a * L^b   (linear in log-log)
    b, loga = np.polyfit(np.log(L), np.log(W), 1)
    a = math.exp(loga)
    pred = a * L**b
    ss_res = float(np.sum((np.log(W) - np.log(pred))**2))
    ss_tot = float(np.sum((np.log(W) - np.log(W).mean())**2))
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")

    # per-structure disk, split into length-scaling cif and ~constant json
    struct_b = dir_bytes(os.path.join(args.test_dir, "structures"))
    metrics_b = dir_bytes(os.path.join(args.test_dir, "metrics"))
    n = len(rows)
    cif_per_res = struct_b / L.sum()          # cif grows with residues
    json_per_seq = metrics_b / n              # metrics json ~constant

    # --- histogram: integrate the fit over the real length distribution ---
    foldable_seqs = foldable_res = skipped_seqs = skipped_res = 0
    gpu_seconds = 0.0
    with open(args.length_hist) as fh:
        for line in fh:
            length, count = line.split()
            length, count = int(length), int(count)
            if length > args.max_len:
                skipped_seqs += count
                skipped_res += length * count
            else:
                foldable_seqs += count
                foldable_res += length * count
                gpu_seconds += count * a * length**b

    gpu_hours = gpu_seconds / 3600
    storage_gb = (cif_per_res * foldable_res + json_per_seq * foldable_seqs) / 1e9

    print(f"=== power-law fit from {n} folds (len {int(L.min())}-{int(L.max())}aa) ===")
    print(f"wall_s = {a:.4g} * L^{b:.3f}   (R^2={r2:.3f} in log space)")
    print(f"  implied: {a*100**b:.2f}s at 100aa, {a*264**b:.2f}s at 264aa (median), {a*2048**b:.1f}s at 2048aa")
    print(f"disk: {cif_per_res:.1f} B/residue (cif) + {json_per_seq:.0f} B/seq (json)")
    print()
    print(f"=== full input ===")
    print(f"foldable (<= {args.max_len}aa): {foldable_seqs:,} seqs, {foldable_res:,} residues")
    print(f"skipped  (>  {args.max_len}aa): {skipped_seqs:,} seqs ({skipped_seqs/(foldable_seqs+skipped_seqs)*100:.2f}%), handle separately")
    print()
    print(f"=== projection (foldable only) ===")
    print(f"total GPU-time : {gpu_hours:,.0f} GPU-hours  ({gpu_hours/24:,.0f} GPU-days)")
    print(f"total storage  : {storage_gb:,.0f} GB")
    fits = storage_gb <= args.disk_budget_gb
    print(f"fits in {args.disk_budget_gb:.0f} GB? {'YES' if fits else 'NO -> fold-wave then sync_to_s3.sh + delete'}")
    if not fits:
        print(f"                 (~{math.ceil(storage_gb/args.disk_budget_gb)} offload waves)")
    print()
    print(f"=== wall-clock by fleet size (100% GPU efficiency) ===")
    print(f"{'GPUs':>6} {'days':>8}")
    for g in (int(x) for x in args.gpu_scenarios.split(",")):
        print(f"{g:>6} {gpu_hours/g/24:>8.1f}")
    print("(real wall-clock is longer: queue waits, per-job model load ~1-2 min, preemption.)")


if __name__ == "__main__":
    main()
