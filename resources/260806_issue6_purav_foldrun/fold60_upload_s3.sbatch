#!/bin/bash
# Upload the 60% fold to S3 and verify it against the local copy (object counts and bytes per prefix).
# Does NOT delete anything locally. Resumable: rerunning only uploads what is missing or changed.
#   sbatch fold60_upload_s3.sbatch
#SBATCH --job-name=fold60-s3-upload
#SBATCH --account=aip-qmorris
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/fold60-s3-upload-%j.out
set -u -o pipefail
cd "${SLURM_SUBMIT_DIR}"
A="$HOME/aws-cli-v2/bin/aws --profile igem-purav"
DST=s3://petadex-protein-structures/petadex-v1.1-60pid-centroids-structures
echo "=== upload start $(date '+%F %H:%M') on $(hostname)"
$A sts get-caller-identity --query Account --output text || { echo "AWS credentials not valid"; exit 1; }
$A s3 sync data/fold60pid "$DST" --exclude "shard_hours_projection.tsv" --only-show-errors
$A s3 sync data/fold60pid_input "$DST/inputs" --only-show-errors
echo "=== upload done $(date '+%F %H:%M'); verifying"
fail=0
check() {  # local_dir s3_prefix
  local lf lb rf rb
  lf=$(find "$1" -type f | wc -l); lb=$(find "$1" -type f -printf '%s\n' | awk '{s+=$1} END{print s+0}')
  read -r rf rb < <($A s3 ls "$2/" --recursive --summarize | awk '/Total Objects/{o=$3} /Total Size/{s=$3} END{print o+0, s+0}')
  local ok=OK; [[ "$lf" == "$rf" && "$lb" == "$rb" ]] || { ok=MISMATCH; fail=1; }
  printf "%-9s %-45s local %9d files %14d B | s3 %9d objects %14d B\n" "$ok" "$2" "$lf" "$lb" "$rf" "$rb"
}
for t in short_full mid_fast tail_h100; do check "data/fold60pid/$t" "$DST/$t"; done
check data/fold60pid/qc "$DST/qc"
check data/fold60pid_input "$DST/inputs"
$A s3 ls "$DST/README.md" >/dev/null && echo "OK        README.md present" || { echo "MISMATCH  README.md missing"; fail=1; }
[[ $fail == 0 ]] && echo "=== VERIFIED: S3 copy matches local $(date '+%F %H:%M')" || echo "=== VERIFICATION FAILED $(date '+%F %H:%M')"
