#!/bin/bash
# SLURM job-array worker: each array task folds one length-band shard on one L40S.
# Submit via submit_fold.sh (sets INPUT / OUTDIR / NSHARDS and the --array range).
# Conventions match ../../../boltz_petase/main_proteins_sbatch.sh (Killarney / aip-wanglab).
#SBATCH --job-name=fold90pid
#SBATCH --account=aip-wanglab
#SBATCH --nodes=1
#SBATCH --gpus-per-node=l40s:1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00                   # 2 days; cost-balanced shard at NSHARDS=1000 is ~36 GPU-h
#SBATCH --output=logs/fold90pid-%A_%a.out
#SBATCH --array=1-16                        # overridden by submit_fold.sh --array
# Killarney L40S partitions are walltime-banded (b1=3h .. b4=3day .. b5=7day). A 2-day request
# routes to gpubase_l40s_b4. If your site needs it explicit, add: #SBATCH --partition=gpubase_l40s_b4

set -e
set -u

export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# SCRIPTDIR is exported by submit_fold.sh; fall back to the submit directory.
HERE="${SCRIPTDIR:-$SLURM_SUBMIT_DIR}"
REPO="$(cd "$HERE/../.." && pwd)"

# Environment: modules + venv + user-site isolation, then HF cache on /project.
source "$REPO/activate_petabite.sh"
export HF_HOME=/home/purav/projects/aip-qmorris/purav/hf_cache
export PYTHONNOUSERSITE=1

# INPUT / OUTDIR / NSHARDS come from submit_fold.sh (exported into the job env).
: "${INPUT:?set INPUT}"; : "${OUTDIR:?set OUTDIR}"; : "${NSHARDS:?set NSHARDS}"
DTYPE="${DTYPE:-fp32}"   # bf16 errors on this build; fp32 confirmed working
FAST_FLAG=""; [[ "${FAST:-0}" == "1" ]] && FAST_FLAG="--fast"   # ESMFold2-Fast when FAST=1

echo "task ${SLURM_ARRAY_TASK_ID}/${NSHARDS} on $(hostname)  (dtype=$DTYPE fast=${FAST:-0})"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

python -u "$HERE/fold.py" "$INPUT" \
    --outdir "$OUTDIR" \
    --shard "${SLURM_ARRAY_TASK_ID}/${NSHARDS}" \
    --num-sampling-steps 68 --num-loops 10 \
    --max-len 1280 \
    --dtype "$DTYPE" $FAST_FLAG
