#!/usr/bin/env bash
# Launch the folding fleet as a SLURM job array. Set the fleet size here after the
# 100-seq test + estimate_fleet.py tell you how many GPUs / how long you need.
#
#   ./submit_fold.sh
#
# Each array task folds shard TASK_ID/NSHARDS. Re-running is safe: finished structures
# (structures/<id>.cif) are skipped, so a partially-done shard resumes where it stopped.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- fleet configuration ----
export INPUT="$HERE/data/90pid_catalytic_orfs.orfid.fa"
export OUTDIR="$HERE/data/fold90pid"      # on /project (740 GB); holds structures/ + metrics/
export NSHARDS=1000                        # max array size on Killarney is 1001; balanced shards ~36 GPU-h each
export DTYPE=fp32                          # bf16 errors on this build; fp32 confirmed working
export FAST=1                              # 1 = ESMFold2-Fast (1.78x faster, -2 pLDDT); 0 = full
MAXCONCURRENT=16                           # simultaneous running tasks; raise to your allocation's GPU limit
export SCRIPTDIR="$HERE"
# -----------------------------

mkdir -p "$OUTDIR" "$HERE/logs"
cd "$HERE"   # so the sbatch's logs/ output path and SLURM_SUBMIT_DIR resolve here

# TESTN=2 ./submit_fold.sh  -> submit only shards 1..TESTN first, to validate the plumbing
ARRAY="1-${NSHARDS}%${MAXCONCURRENT}"
[[ -n "${TESTN:-}" ]] && ARRAY="1-${TESTN}"

echo "submitting array ${ARRAY}  INPUT=$INPUT  OUTDIR=$OUTDIR  FAST=$FAST"
sbatch --array="$ARRAY" \
       --export=ALL,INPUT,OUTDIR,NSHARDS,DTYPE,FAST,SCRIPTDIR \
       "$HERE/fold_shard.sbatch"
