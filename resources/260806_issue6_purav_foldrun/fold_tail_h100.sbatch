#!/bin/bash
# Long tail (> 1740 aa, data/fold60pid_input/tail_h100.fa, 3,440 for the 60% level) on one H100.
# NOT YET TESTED: sequences beyond what an 80 GB H100 fits even with chunk 32 will OOM and be
# logged; those (>~4000 aa) need Fold-CP across GPUs or another strategy. Run a TESTN-sized
# trial first:  sbatch --array=1-1 fold_tail_h100.sbatch
#SBATCH --job-name=fold60-tail
#SBATCH --account=aip-qmorris
#SBATCH --nodes=1
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/fold60-tail-%A_%a.out
set -u -o pipefail
export PYTHONUNBUFFERED=1 PYTHONFAULTHANDLER=1 HF_HUB_OFFLINE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd "${SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1
python fold.py data/fold60pid_input/tail_h100.fa --outdir data/fold60pid/tail_h100 \
    --shard "${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}" \
    --backend native --fast --esmc-precision bf16 --kernel fused --chunk 32 \
    --num-sampling-steps 32 --num-loops 3 --seed 0 --max-len 100000 \
    --cost-exp 1.66 --cost-floor 1000 --gzip --subdirs 1000 \
    2>&1 | grep -vE "Warning|warn\("
