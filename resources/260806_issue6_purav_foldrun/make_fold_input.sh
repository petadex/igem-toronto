#!/usr/bin/env bash
# Build the ORFid-keyed FASTA that fold.py consumes.
#
# Rewrites the catalytic-core FASTA headers so each record is keyed by its ORFid:
#   >3|7QJM|||||/95-371   ->   >orf3
# fold.py uses that >orf<id> as the output filename, so structures join back on ORFid.
#
# Provenance of the input 90pid_catalytic_orfs.fa (already on disk), for reproducibility:
#   1. ID list (18,172,960 rows) from the DB:
#        SELECT DISTINCT "centroid_orf_id" FROM "90pid_family_clusters";
#   2. Pulled the matching sequences out of the 108 GiB catalytic-ORF FASTA:
#        aws s3 cp s3://petadex/logan/petadex.catalytic_orfs.v1.1.fa - --no-sign-request \
#          | seqkit grep --id-regexp '^([^|]+)' -f 90_orfs.txt > 90pid_catalytic_orfs.fa
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
IN="${1:-$REPO/90pid_catalytic_orfs.fa}"
OUTDIR="$HERE/data"
OUT="$OUTDIR/90pid_catalytic_orfs.orfid.fa"

mkdir -p "$OUTDIR"

if [[ -s "$OUT" ]]; then
    echo "already built: $OUT (delete to rebuild)"
else
    echo "rewriting headers: $IN -> $OUT"
    awk '/^>/{split(substr($0,2),a,"|"); print ">orf" a[1]; next} {print}' "$IN" > "$OUT"
fi

echo "=== seqkit stats (num_seqs should be 18,172,960; sum_len = total residues) ==="
seqkit stats -T "$OUT"
