#!/bin/bash
# Phase 0: software-stack benchmark. Native esm 3.4.1 EsmFold2Model vs the Aug 2026 HF-fork
# folds (data/sample = Full, data/sample_fast = Fast, both fp32, 68 steps / 10 loops) on the
# identical 79-sequence length-stratified sample (--sample 80 --max-len 2048, seed 0).
# All steps skip finished .cif files, so resubmitting resumes.
#SBATCH --job-name=fold-stack
#SBATCH --account=aip-qmorris
#SBATCH --nodes=1
#SBATCH --gpus-per-node=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=03:00:00
#SBATCH --output=logs/stack-%j.out
set -u
export PYTHONUNBUFFERED=1 PYTHONFAULTHANDLER=1 HF_HUB_OFFLINE=1
cd "${SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

IN=data/90pid_catalytic_orfs.orfid.fa
S="--sample 80 --max-len 2048 --backend native"
run() { echo "=== $1"; shift; python fold.py "$IN" $S "$@" 2>&1 | grep -E "loaded|complete|FAILED|Error"; }

run "fast bf16 kernel=none"           --outdir data/stack/fast_bf16_none  --fast --kernel none           --chunk none
run "fast bf16 kernel=fused"          --outdir data/stack/fast_bf16_fused --fast --kernel fused          --chunk none
run "fast bf16 kernel=cueq"           --outdir data/stack/fast_bf16_cueq  --fast --kernel cuequivariance --chunk none
run "fast bf16 fused seed=1 (noise)"  --outdir data/stack/fast_bf16_fused_seed1 --fast --kernel fused --chunk none --seed 1
run "fast esmc-fp32 fused"            --outdir data/stack/fast_esmcfp32_fused --fast --kernel fused --chunk none --esmc-precision fp32
run "full bf16 fused"                 --outdir data/stack/full_bf16_fused --kernel fused --chunk none
run "full bf16 fused seed=1 (noise)"  --outdir data/stack/full_bf16_fused_seed1 --kernel fused --chunk none --seed 1

# memory ceiling with chunking, both models (ladder 505-2052aa)
echo "=== ladder fast chunk64";  python fold.py data/probe_ladder.fa --outdir data/stack/ladder_fast_c64 --backend native --fast --kernel fused --chunk 64 --max-len 100000 2>&1 | grep -E "len=|FAILED"
echo "=== ladder full chunk64";  python fold.py data/probe_ladder.fa --outdir data/stack/ladder_full_c64 --backend native --kernel fused --chunk 64 --max-len 100000 2>&1 | grep -E "len=|FAILED"
echo "=== DONE"
