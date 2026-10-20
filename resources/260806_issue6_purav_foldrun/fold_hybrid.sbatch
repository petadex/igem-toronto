#!/bin/bash
# Production hybrid fold (Part II, Step 15): one array task folds one cost-balanced shard of one tier.
#   short_full : ESMFold2 Full,  32 steps / 3 loops, L < 150
#   mid_fast   : ESMFold2-Fast,  32 steps / 3 loops, 150 <= L <= 1740
# Native esm 3.4.1, ESMC bf16, fused kernels, auto chunking, gzip mmCIF in 1000 subdirectories.
# 3 h walltime -> gpubase_l40s_b1 (largest L40S pool). A task cut off at 3 h resumes on resubmit
# (finished structures are skipped). Submit with submit_fold60.sh, which sets TIER/INPUT/OUTDIR/NSHARDS.
#SBATCH --job-name=fold60
#SBATCH --account=aip-qmorris
#SBATCH --nodes=1
#SBATCH --gpus-per-node=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=03:00:00
#SBATCH --output=logs/fold60-%x-%A_%a.out
set -u -o pipefail
export PYTHONUNBUFFERED=1 PYTHONFAULTHANDLER=1 HF_HUB_OFFLINE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd "${SCRIPTDIR:-$SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1

: "${TIER:?}"; : "${INPUT:?}"; : "${OUTDIR:?}"; : "${NSHARDS:?}"
case "$TIER" in
  short_full) MODEL="";       COST="--cost-exp 0.16 --cost-floor 0";    MAXLEN=149 ;;
  mid_fast)   MODEL="--fast"; COST="--cost-exp 1.66 --cost-floor 1000"; MAXLEN=1740 ;;
  *) echo "unknown TIER=$TIER"; exit 2 ;;
esac
echo "tier=$TIER shard ${SLURM_ARRAY_TASK_ID}/${NSHARDS} on $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# Preflight. A node whose GPU is visible to nvidia-smi but not to CUDA (kn140, 2026-10-15) makes
# every task fail in ~16 s, and the instantly freed slot pulls the next task onto the same node
# (43 shards lost there in 20 min). Hold the slot for 15 min before failing, so the array drains
# slowly, and print the node so fold60_watch.sh / a human can exclude it.
if ! python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null; then
  echo "PREFLIGHT FAILED: CUDA unavailable on $(hostname); holding 15 min, then exiting 1"
  sleep 900
  exit 1
fi

python fold.py "$INPUT" --outdir "$OUTDIR" --shard "${SLURM_ARRAY_TASK_ID}/${NSHARDS}" \
    --backend native $MODEL --esmc-precision bf16 --kernel fused --chunk auto \
    --num-sampling-steps 32 --num-loops 3 --seed 0 --max-len "$MAXLEN" \
    $COST --gzip --subdirs 1000 \
    2>&1 | grep -vE "Warning|warn\("
