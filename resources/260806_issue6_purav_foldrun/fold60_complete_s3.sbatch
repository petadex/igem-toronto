#!/bin/bash
# Complete the 60% structure set on S3 (see fold60_complete.py):
#   - copy the 19,430 AFDB structures into afdb/structures/<NNN>/orf<ID>.cif.gz (sequence-checked, gzipped)
#   - fold-instead: 27 AFDB keys are not the exact sequence (qc/afdb_inexact_refolded.tsv); their ESMFold2
#     structures (folded into data/fold60pid/mid_fast by fold_hybrid.sbatch) replace the AFDB copies
#   - copy the 31 duplicates' twin structure + metrics under the duplicate's ID (server-side)
#   - upload structure_index.tsv, the updated README and the QC lists; verify everything is there
# Deletes on S3 only the 27 inexact AFDB copies (listed in qc/afdb_remove.tsv). Resumable.
#   sbatch fold60_complete_s3.sbatch
#SBATCH --job-name=fold60-s3-complete
#SBATCH --account=aip-qmorris
#SBATCH --cpus-per-task=32
#SBATCH --mem=32gb
#SBATCH --time=6:00:00
#SBATCH --output=logs/fold60-s3-complete-%j.out
set -u -o pipefail
cd "${SLURM_SUBMIT_DIR}"
# one SSO lookup for the whole job: each `aws --profile` call otherwise re-fetches role credentials (16 s vs 5 s per file)
eval "$($HOME/aws-cli-v2/bin/aws configure export-credentials --profile igem-purav --format env)" || { echo "AWS credentials not valid"; exit 1; }
export AWS_DEFAULT_REGION=us-east-1  # bucket region
export AWS="$HOME/aws-cli-v2/bin/aws"
DST=s3://petadex-protein-structures/petadex-v1.1-60pid-centroids-structures
Q=data/fold60pid/qc
export S="/scratch/$USER/fold60_afdb"  # on /scratch so a rerun resumes
mkdir -p "$S" "$Q"
echo "=== start $(date '+%F %H:%M') on $(hostname)"
$AWS sts get-caller-identity --query Account --output text || { echo "AWS credentials not valid"; exit 1; }

python3 fold60_complete.py lists "$Q" || exit 1
mv "$Q/structure_index.tsv" data/fold60pid/structure_index.tsv
n=$(grep -c '^>' data/fold60pid_input/afdb_refold.fa)
k=$(grep '^>' data/fold60pid_input/afdb_refold.fa | sed 's/^>orf//' | while read -r o; do
      ls data/fold60pid/mid_fast/structures/$(printf %03d $((o % 1000)))/orf$o.cif.gz 2>/dev/null; done | wc -l)
echo "refolded structures present: $k of $n"; [ "$k" -eq "$n" ] || { echo "refold incomplete; stopping"; exit 1; }
# drop the inexact AFDB files from staging so the sync below cannot upload them again
sed "s|^$DST/||" "$Q/afdb_remove.tsv" | while read -r r; do rm -f "$S/$r" "$S/${r%.gz}"; done

echo "=== AFDB download $(date '+%H:%M')"
fetch() { [ -s "$S/$2" ] || { mkdir -p "$S/$(dirname "$2")" && $AWS s3 cp "$1" "$S/$2" --only-show-errors; }; }
export -f fetch
for try in 1 2 3; do  # second pass retries anything that failed
  cut -f1,2 "$Q/afdb_jobs.tsv" | xargs -P 64 -n 2 bash -c 'fetch "$1" "$2"' _
done
python3 fold60_complete.py check-afdb "$S" "$Q" || { echo "AFDB files missing; stopping"; exit 1; }

echo "=== gzip + upload AFDB $(date '+%H:%M')"
find "$S/afdb" -name '*.cif' -print0 | xargs -0 -P 32 -n 50 gzip -n -f
$AWS s3 sync "$S/afdb" "$DST/afdb" --only-show-errors
echo "=== replace the 27 inexact AFDB copies with the ESMFold2 folds $(date '+%H:%M')"
$AWS s3 sync data/fold60pid/mid_fast "$DST/mid_fast" --only-show-errors
while read -r k; do $AWS s3 rm "$k" --only-show-errors; done < "$Q/afdb_remove.tsv"

echo "=== duplicates $(date '+%H:%M')"
xargs -P 16 -n 2 $AWS s3 cp --only-show-errors < "$Q/dup_jobs.tsv"

echo "=== index, README, qc $(date '+%H:%M')"
$AWS s3 cp data/fold60pid/structure_index.tsv "$DST/structure_index.tsv" --only-show-errors
$AWS s3 cp data/fold60pid/README.md "$DST/README.md" --only-show-errors
$AWS s3 sync "$Q" "$DST/qc" --only-show-errors

echo "=== verify $(date '+%H:%M')"
fail=0
want_afdb=$(( $(wc -l < "$Q/afdb_jobs.tsv") + $(grep -c "/afdb/" "$Q/dup_jobs.tsv") ))
n=$($AWS s3 ls "$DST/afdb/" --recursive | grep -c '\.cif\.gz$')
echo "afdb objects: $n (expect $want_afdb)"; [ "$n" -eq "$want_afdb" ] || fail=1
r=0; while read -r k; do $AWS s3 ls "$k" > /dev/null && r=$((r + 1)); done < "$Q/afdb_remove.tsv"
echo "inexact AFDB copies still on S3: $r (expect 0)"; [ "$r" -eq 0 ] || fail=1
m=0
while read -r src dst; do $AWS s3 ls "$dst" > /dev/null && m=$((m + 1)); done < "$Q/dup_jobs.tsv"
echo "duplicate objects: $m (expect 61)"; [ "$m" -eq 61 ] || fail=1
# every structure_key in the index exists on S3
$AWS s3 ls "$DST/" --recursive | awk '{print $4}' | sed "s|^petadex-v1.1-60pid-centroids-structures/||" | sort > "$S/keys.txt"
awk -F'\t' 'NR>1 && $3!=""{print $3} NR>1 && $4!=""{print $4}' data/fold60pid/structure_index.tsv | sort > "$S/want.txt"
miss=$(comm -13 "$S/keys.txt" "$S/want.txt" | wc -l)
echo "index entries missing on S3: $miss of $(wc -l < "$S/want.txt")"; [ "$miss" -eq 0 ] || fail=1
[ $fail -eq 0 ] && echo "=== VERIFIED: S3 set complete $(date '+%F %H:%M')" || echo "=== VERIFY FAILED"
exit $fail
