#!/bin/bash
set -euo pipefail

IPRSCAN=~/interproscan/interproscan-5.72-103.0/interproscan.sh
CHUNKS_DIR=~/iprscan_chunks
RESULTS_DIR=~/iprscan_results
APPL="Pfam,TIGRFAM,CDD,SUPERFAMILY,Gene3D,PANTHER"
CPU=8

mkdir -p "$RESULTS_DIR"

START_CHUNK=1   # first chunk to process (inclusive)
END_CHUNK=13    # last chunk to process (inclusive)

total=$(ls "$CHUNKS_DIR"/*.fasta | wc -l)
count=0

for f in "$CHUNKS_DIR"/*.fasta; do
    count=$((count + 1))
    base=$(basename "$f" .fasta)
    out="$RESULTS_DIR/$base"

    # Skip chunks outside the requested range
    if [ "$count" -lt "$START_CHUNK" ] || [ "$count" -gt "$END_CHUNK" ]; then
        echo "[$count/$total] Skipping $base (outside range $START_CHUNK-$END_CHUNK)"
        continue
    fi

    # Skip if already done (with or without .tsv extension)
    if [ -f "${out}.tsv" ] || [ -f "${out}" ]; then
        echo "[$count/$total] Skipping $base (already done)"
        continue
    fi

    echo "[$count/$total] Running $base ..."
    "$IPRSCAN" \
        -i "$f" \
        -o "$out" \
        -f TSV \
        -goterms \
        -appl "$APPL" \
        -cpu "$CPU" \
        -dp \
        --disable-precalc 2>&1 | tail -5

    echo "[$count/$total] Done: $base"
done

echo ""
echo "All chunks complete. Merging TSV results..."
# Gather all result files (with or without .tsv extension), excluding the merged output
mapfile -t result_files < <(ls "$RESULTS_DIR"/chunk_*.tsv "$RESULTS_DIR"/chunk_* 2>/dev/null | grep -v "all_results" | sort -u)
cat "${result_files[@]}" > "$RESULTS_DIR/all_results.tsv"
echo "Merged output: $RESULTS_DIR/all_results.tsv"
echo "Total lines: $(wc -l < $RESULTS_DIR/all_results.tsv)"
