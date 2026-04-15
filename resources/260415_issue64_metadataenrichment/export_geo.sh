#!/bin/bash

CONN="postgresql://public_reader:serratus@serratus-aurora-20210406.cluster-ro-ccz9y6yshbls.us-east-1.rds.amazonaws.com:5432/logan"

QUERY="SELECT DISTINCT accession, ST_Y(lat_lon) AS latitude, ST_X(lat_lon) AS longitude, elevation, country, biome, confidence FROM biosample_geographical_location"

export_chunk() {
    local label=$1
    local where_clause=$2
    local outfile=$3

    if [ -f "$outfile" ]; then
        echo "Skipping $label — $outfile already exists"
        return
    fi

    echo "Exporting $label..."
    psql "$CONN" -c "\copy ($QUERY WHERE $where_clause) TO '$outfile' CSV HEADER"

    if [ $? -eq 0 ]; then
        echo "  Done: $(wc -l < "$outfile") rows"
    else
        echo "  FAILED: $label"
    fi
}

# SAMD
export_chunk "SAMD" "accession LIKE 'SAMD%'" "geo_samd.csv"

# SAMEA 2-9
for i in 2 3 4 5 6 7 8 9; do
    export_chunk "SAMEA${i}" "accession LIKE 'SAMEA${i}%'" "geo_samea${i}.csv"
done

# SAMEA1 - split further
for i in 0 2 3 4; do
    export_chunk "SAMEA1${i}" "accession LIKE 'SAMEA1${i}%'" "geo_samea1${i}.csv"
done

# SAMEA15-19 combined (small)
export_chunk "SAMEA15-19" "accession >= 'SAMEA15' AND accession < 'SAMEA2'" "geo_samea15-19.csv"

# SAMEA11 - deeper split
for i in 0 1 2 3 4 5 6 7 8 9; do
    export_chunk "SAMEA11${i}" "accession LIKE 'SAMEA11${i}%'" "geo_samea11${i}.csv"
done

# SAMN00-41
for i in $(seq -w 0 41); do
    export_chunk "SAMN${i}" "accession LIKE 'SAMN${i}%'" "geo_samn${i}.csv"
done

echo "All done!"