#!/usr/bin/env bash
# Probe public readability of Alex structure objects (Range GET).
# Usage: ./probe_structure_s3.sh [orfId]
set -euo pipefail

ORF_ID="${1:-4981589}"
BASE="https://petadex-protein-structures.s3.amazonaws.com"

paths=(
  "esmfold2-centroids/test2/structures/orf${ORF_ID}.cif"
  "esmfold2-centroids/test2/metrics/orf${ORF_ID}.json"
  "esmfold2-centroids/60pid/structures/orf${ORF_ID}.cif"
  "esmfold2-centroids/60pid/metrics/orf${ORF_ID}.json"
  "esmfold2-centroids/60pid-msa/structures/orf${ORF_ID}.cif"
  "esmfold2-centroids/60pid-msa/metrics/orf${ORF_ID}.json"
)

echo "Probing ORF ${ORF_ID}"
echo

for path in "${paths[@]}"; do
  code=$(curl -sS -o /dev/null -w "%{http_code}" -H "Range: bytes=0-0" "${BASE}/${path}" || echo "ERR")
  printf "%s  %s\n" "$code" "$path"
done

echo
echo "206/200 = public. 403 = ACL closed. 404 = missing key."
