#!/bin/bash
# Phase 0 smoke: native esm 3.4.1 ESMFold2-Fast loads from the cached weights and folds
# 3 x ~260aa with each kernel backend. Run inside a GPU allocation.
set -u
cd "$(dirname "$0")"
source ../../activate_esmfold2.sh >/dev/null 2>&1
export HF_HUB_OFFLINE=1 PYTHONUNBUFFERED=1
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
for K in none fused cuequivariance; do
  echo "=== kernel=$K"
  python fold.py data/smoke3.fa --outdir data/smoke_native_$K --backend native --fast \
      --esmc-precision bf16 --kernel $K --chunk none 2>&1 | grep -vE "Warning|warn\(" | tail -6
done
