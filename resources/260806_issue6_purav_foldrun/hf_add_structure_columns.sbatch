#!/bin/bash
# Part III: rewrite all 131 orfs/ parts of petadex/petadex-corpus-full with is_c60_centroid + structure_path,
# check every part and the totals, then open a Hugging Face PR (never pushes to main). Resumable.
#   sbatch hf_add_structure_columns.sbatch
#SBATCH --job-name=hf-corpus-structcols
#SBATCH --account=aip-qmorris
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=8:00:00
#SBATCH --output=logs/hf-corpus-structcols-%j.out
set -eu -o pipefail
cd "${SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1
echo "=== start $(date '+%F %H:%M') on $(hostname)"
python hf_add_structure_columns.py lookup
python hf_add_structure_columns.py rewrite
echo "=== rewrite done $(date '+%H:%M')"
python hf_add_structure_columns.py totals
python hf_add_structure_columns.py pr
echo "=== done $(date '+%F %H:%M')"
