#!/bin/bash
# Phase 2: model x diffusion-steps x refinement-loops sweep on benchmark set B
# (data/bench/bench_B.fa, 620 stratified 90pid centroids). One array task per row of
# sweep_configs.tsv. Native esm 3.4.1 stack, ESMC bf16, fused kernels, auto chunking.
# Reference = Full, 100 steps, 20 loops. Seeds 1/2 on the 68/10 configs give the noise floor.
#   sbatch --array=1-17%8 sweep.sbatch                                     # set B
#   INPUT=data/gt/gt.fa OUTROOT=data/sweep_gt sbatch --array=1-17%8 sweep.sbatch   # PDB ground truth
#SBATCH --job-name=fold-sweep
#SBATCH --account=aip-qmorris
#SBATCH --nodes=1
#SBATCH --gpus-per-node=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/sweep-%A_%a.out
set -u -o pipefail   # a python failure must fail the task, not be masked by the grep pipe
export PYTHONUNBUFFERED=1 PYTHONFAULTHANDLER=1 HF_HUB_OFFLINE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd "${SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1

read -r _ NAME MODEL STEPS LOOPS SEED < <(awk -F'\t' -v t="$SLURM_ARRAY_TASK_ID" '$1==t' sweep_configs.tsv)
: "${NAME:?no config row for task $SLURM_ARRAY_TASK_ID}"
FAST=""; [[ "$MODEL" == "fast" ]] && FAST="--fast"
echo "task $SLURM_ARRAY_TASK_ID: $NAME  model=$MODEL steps=$STEPS loops=$LOOPS seed=$SEED  $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

python fold.py "${INPUT:-data/bench/bench_B.fa}" --outdir "${OUTROOT:-data/sweep}/$NAME" --backend native $FAST \
    --esmc-precision bf16 --kernel fused --chunk auto --chunk-above "${CHUNK_ABOVE:-1000}" \
    --num-sampling-steps "$STEPS" --num-loops "$LOOPS" --seed "$SEED" --max-len 2048 \
    2>&1 | grep -vE "Warning|warn\("
