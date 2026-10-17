#!/bin/bash
# Delete the local 60% structures and metrics after the verified S3 upload (fold60_upload_s3.sbatch, job 6323482).
# Keeps progress CSVs, qc/, README.md and data/fold60pid_input/. Run only after the user has confirmed.
#   sbatch fold60_delete_local.sbatch
#SBATCH --job-name=fold60-delete-local
#SBATCH --account=aip-qmorris
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=6:00:00
#SBATCH --output=logs/fold60-delete-local-%j.out
set -u -o pipefail
cd "${SLURM_SUBMIT_DIR}"
echo "=== delete start $(date '+%F %H:%M') on $(hostname)"
for t in short_full mid_fast tail_h100; do
  for d in structures metrics; do
    p="data/fold60pid/$t/$d"
    [ -d "$p" ] || { echo "$p: already gone"; continue; }
    # one rm per <NNN> subdirectory, 8 in parallel
    find "$p" -mindepth 1 -maxdepth 1 -type d -print0 | xargs -0 -P 8 -n 20 rm -rf
    rmdir "$p" && echo "$p: deleted $(date '+%H:%M')"
  done
done
echo "=== delete done $(date '+%F %H:%M')"
du -sh data/fold60pid data/fold60pid_input
find data/fold60pid | wc -l
