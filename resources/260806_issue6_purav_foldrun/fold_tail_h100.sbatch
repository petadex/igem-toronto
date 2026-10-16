#!/bin/bash
# Long tail on H100: 1741-2500 aa (above the L40S ceiling for Fast, at or below the 2500 aa hard
# cap). ESMFold2-Fast 32 steps / 3 loops like the >=150 aa tier, chunk 32 so the L^2 pair ops fit
# 80 GB at 2500 aa. Sequences > 2500 aa are not folded (hard cap, decision 2026-10-16).
# Shard NSHARDS holds the longest sequences, so run it alone first as the memory test:
#   NSHARDS=6 sbatch --array=6 --export=ALL,NSHARDS=6 fold_tail_h100.sbatch
#   NSHARDS=6 sbatch --array=1-5 --export=ALL,NSHARDS=6 fold_tail_h100.sbatch
#SBATCH --job-name=fold60-tail
#SBATCH --account=aip-qmorris
#SBATCH --nodes=1
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=03:00:00
#SBATCH --output=logs/fold60-tail-%A_%a.out
set -u -o pipefail
export PYTHONUNBUFFERED=1 PYTHONFAULTHANDLER=1 HF_HUB_OFFLINE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd "${SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1
: "${NSHARDS:?}"
echo "tail shard ${SLURM_ARRAY_TASK_ID}/${NSHARDS} on $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
if ! python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null; then
  echo "PREFLIGHT FAILED: CUDA unavailable on $(hostname); holding 15 min, then exiting 1"; sleep 900; exit 1
fi
python fold.py data/fold60pid_input/tail_h100_le2500.fa --outdir data/fold60pid/tail_h100 \
    --shard "${SLURM_ARRAY_TASK_ID}/${NSHARDS}" \
    --backend native --fast --esmc-precision bf16 --kernel fused --chunk 32 \
    --num-sampling-steps 32 --num-loops 3 --seed 0 --max-len 2500 \
    --cost-exp 1.66 --cost-floor 1000 --gzip --subdirs 1000 \
    2>&1 | grep -vE "Warning|warn\("
