#!/bin/bash
# Split a FASTA file into chunks of N sequences each
INPUT=$1
OUTDIR=$2
CHUNK_SIZE=${3:-5000}

mkdir -p "$OUTDIR"

awk -v outdir="$OUTDIR" -v chunk="$CHUNK_SIZE" '
/^>/ {
    count++
    if ((count-1) % chunk == 0) {
        filenum = int((count-1) / chunk) + 1
        close(outfile)
        outfile = outdir "/chunk_" sprintf("%04d", filenum) ".fasta"
    }
    print > outfile
    next
}
{ print > outfile }
' "$INPUT"

echo "Done. $(ls $OUTDIR/*.fasta | wc -l) chunks written to $OUTDIR"
