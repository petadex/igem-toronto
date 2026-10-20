#!/bin/bash
# One-shot calibration job: everything needed to choose full-vs-Fast and the fp32 memory
# ceiling, run unattended on one L40S. Submit from this dir:  sbatch calibrate.sbatch
# All steps are skip-if-done (finished .cif files are skipped), so it is safe to resubmit.
#SBATCH --job-name=fold-calibrate
#SBATCH --account=aip-wanglab
#SBATCH --nodes=1
#SBATCH --gpus-per-node=l40s:1
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/calibrate-%j.out

set -e
set -u
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

HERE="${SLURM_SUBMIT_DIR}"
cd "$HERE"
REPO="$(cd "$HERE/../.." && pwd)"
source "$REPO/activate_petabite.sh"
export HF_HOME=/home/purav/projects/aip-qmorris/purav/hf_cache
export PYTHONNOUSERSITE=1

IN=data/90pid_catalytic_orfs.orfid.fa
HIST=data/length_hist.tsv

nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# 0) length-ladder probe fasta: up to 2 sequences at each target length, to find the
#    fp32 OOM ceiling and get high-L timing for BOTH models on identical inputs.
LADDER=data/probe_ladder.fa
if [[ ! -s "$LADDER" ]]; then
  : > "$LADDER"
  for L in 512 768 1024 1280 1536 1792 2048; do
    seqkit seq -m $((L-8)) -M $((L+8)) "$IN" | seqkit head -n 2 >> "$LADDER" || true
  done
fi
echo "=== ladder probe sequences ==="; seqkit fx2tab -nil "$LADDER" || true

# 1) full-model stratified sample (already done -> skipped)
echo "=== [full] stratified sample (data/sample) ==="
python fold.py "$IN" --outdir data/sample --sample 80 --max-len 2048

# 2) Fast-model stratified sample (identical 80 seqs, seed 0)
echo "=== [fast] stratified sample (data/sample_fast) ==="
python fold.py "$IN" --outdir data/sample_fast --sample 80 --max-len 2048 --fast

# 3) memory-ceiling + high-L timing, both models. max-len huge so nothing is skipped;
#    per-sequence OOMs are caught by fold.py and logged, so the job continues. The first
#    length that errors marks each model's fp32 ceiling.
echo "=== [full] ladder probe (data/probe_full) ==="
python fold.py "$LADDER" --outdir data/probe_full --max-len 100000 || true
echo "=== [fast] ladder probe (data/probe_fast) ==="
python fold.py "$LADDER" --outdir data/probe_fast --max-len 100000 --fast || true

# 4) fleet estimates for both models
echo "=== [full] estimate ==="
python estimate_fleet.py --test-dir data/sample      --length-hist "$HIST" --max-len 2048 || true
echo "=== [fast] estimate ==="
python estimate_fleet.py --test-dir data/sample_fast --length-hist "$HIST" --max-len 2048 || true

# 5) quality + speed comparison on the identical stratified sequences
echo "=== full vs Fast: quality + speed on shared seqs ==="
python - <<'PY' || true
import json, glob, statistics as st
def load(d):
    m = {}
    for f in glob.glob(f"{d}/metrics/*.json"):
        r = json.load(open(f)); m[r["id"]] = r
    return m
full, fast = load("data/sample"), load("data/sample_fast")
ids = sorted(set(full) & set(fast))
if ids:
    print(f"shared sequences: {len(ids)}")
    mp_f = st.mean(full[i]["mean_plddt"] for i in ids)
    mp_s = st.mean(fast[i]["mean_plddt"] for i in ids)
    pt_f = st.mean(full[i]["ptm"] for i in ids)
    pt_s = st.mean(fast[i]["ptm"] for i in ids)
    w_f = st.median(full[i]["wall_s"] for i in ids)
    w_s = st.median(fast[i]["wall_s"] for i in ids)
    print(f"mean pLDDT  full={mp_f:.1f}  fast={mp_s:.1f}  (delta {mp_f-mp_s:+.2f})")
    print(f"mean pTM    full={pt_f:.3f} fast={pt_s:.3f} (delta {pt_f-pt_s:+.3f})")
    print(f"median wall full={w_f:.2f}s fast={w_s:.2f}s  (fast speedup {w_f/w_s:.2f}x)")
else:
    print("no shared sequences yet")
PY

echo "=== DONE calibrate ==="
