#!/usr/bin/env bash
#
# run_clusters.sh -- run the full sequence-purchasing pipeline over a folder of
# aligned-core FASTAs, one cluster per file.  Preprocessing is NOT included:
# every input must already be an aligned core FASTA (gaps as '-', headers
# >coreN_n<k>), exactly like stage0_cluster1/.../c1.core.aln.fasta.
#
# For each <name>.fasta it runs, into <out>/<name>/:
#     1. cutsearch_design.py     the design (writes a timestamped run dir)
#     2. validate_design.py      18 independent checks
#     3. annotate_degenerate.py  the degenerate-codon annotation
#     4. make_order.py           order-ready constructs + digest/ligate round trip
#     5. make_final_figures.py   the result figures        (only with --figures)
#
# Each cluster's output is teed to <out>/<name>/<name>.log, and a one-line-per-
# cluster summary is printed and written to <out>/summary.tsv.  A cluster that
# fails does not stop the batch.
#
# USAGE
#   ./run_clusters.sh <fasta-dir> [out-dir] [--figures]
#
# ENV OVERRIDES (design parameters default to the cluster-1 v2 run)
#   PYTHON        interpreter (default: python3; needs 3.9+, and matplotlib only
#                 for --figures).  To use this machine's conda env from WSL:
#                 PYTHON=/mnt/c/Users/27leo/miniconda3/envs/igem-alex/python.exe
#   FASTA_GLOB    filename pattern within <fasta-dir> (default: *.fasta)
#   MAX_LIBRARY   --max-library      (default 2500)
#   MAX_NT        --max-nt           (default 45000)
#   K_MAX         --k-max            (default 10)
#   SEED          --seed             (default 0)
#   MAX_OLIGO_NT  --max-oligo-nt     (default: unset = no fragment-length cap;
#                 set to 330 for the 350 bp ordered-construct limit)
#   CHEMISTRY     --chemistry        (default gg)
#   CUT_SEARCH    --cut-search       (default dp)
#
set -uo pipefail

# ---- arguments ------------------------------------------------------------- #
FASTA_DIR="${1:-}"
OUT_DIR="${2:-./cluster_runs}"
FIGURES=0
for arg in "$@"; do [ "$arg" = "--figures" ] && FIGURES=1; done

if [ -z "$FASTA_DIR" ] || [ ! -d "$FASTA_DIR" ]; then
    echo "usage: $0 <fasta-dir> [out-dir] [--figures]" >&2
    exit 2
fi

# ---- locate the pipeline scripts (independent of the caller's cwd) ---------- #
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FINAL="$HERE/final"

PYTHON="${PYTHON:-python3}"
FASTA_GLOB="${FASTA_GLOB:-*.fasta}"
MAX_LIBRARY="${MAX_LIBRARY:-2500}"
MAX_NT="${MAX_NT:-45000}"
K_MAX="${K_MAX:-10}"
SEED="${SEED:-0}"
CHEMISTRY="${CHEMISTRY:-gg}"
CUT_SEARCH="${CUT_SEARCH:-dp}"

# --max-oligo-nt is only passed when MAX_OLIGO_NT is set, so "unset" means no cap
CAP_ARGS=()
[ -n "${MAX_OLIGO_NT:-}" ] && CAP_ARGS=(--max-oligo-nt "$MAX_OLIGO_NT")

mkdir -p "$OUT_DIR"
SUMMARY="$OUT_DIR/summary.tsv"
printf 'cluster\tstatus\trecommended_K\tcoverage\tnt_ordered\trun_dir\n' > "$SUMMARY"

echo "pipeline scripts : $FINAL"
echo "python           : $PYTHON  ($("$PYTHON" --version 2>&1))"
echo "input FASTAs     : $FASTA_DIR/$FASTA_GLOB"
echo "output           : $OUT_DIR"
echo "design params    : --chemistry $CHEMISTRY --cut-search $CUT_SEARCH" \
     "--max-library $MAX_LIBRARY --max-nt $MAX_NT --k-max $K_MAX --seed $SEED" \
     "${CAP_ARGS[*]:-(no length cap)}"
echo "figures          : $([ "$FIGURES" = 1 ] && echo yes || echo no)"
echo

shopt -s nullglob
FASTAS=("$FASTA_DIR"/$FASTA_GLOB)
if [ ${#FASTAS[@]} -eq 0 ]; then
    echo "no files matching $FASTA_GLOB in $FASTA_DIR" >&2
    exit 1
fi

# ---- per-cluster pipeline (streams live; a failure is recorded, not fatal) -- #
for aln in "${FASTAS[@]}"; do
    name="$(basename "$aln")"; name="${name%%.*}"        # strip all extensions
    sub="$OUT_DIR/$name"
    log="$sub/$name.log"
    mkdir -p "$sub"; : > "$log"
    echo "=================================================================="
    echo "=== $name"
    echo "=================================================================="

    ok=1
    "$PYTHON" "$FINAL/cutsearch_design.py" "$aln" \
        --chemistry "$CHEMISTRY" --cut-search "$CUT_SEARCH" \
        --max-library "$MAX_LIBRARY" --max-nt "$MAX_NT" \
        --k-max "$K_MAX" --seed "$SEED" "${CAP_ARGS[@]}" \
        --out-dir "$sub" 2>&1 | tee -a "$log" || ok=0

    run="$(ls -dt "$sub"/*/ 2>/dev/null | head -1)"; run="${run%/}"

    if [ "$ok" = 1 ] && [ -n "$run" ] && [ -f "$run/summary.json" ]; then
        "$PYTHON" "$FINAL/validate_design.py"     "$run" --aln "$aln" 2>&1 | tee -a "$log" || ok=0
        "$PYTHON" "$FINAL/annotate_degenerate.py" "$run"             2>&1 | tee -a "$log" || ok=0
        "$PYTHON" "$FINAL/make_order.py"          "$run"             2>&1 | tee -a "$log" || ok=0
        if [ "$FIGURES" = 1 ]; then
            "$PYTHON" "$HERE/make_final_figures.py" "$run" --cluster "$name" 2>&1 | tee -a "$log" || ok=0
        fi
    else
        ok=0
    fi

    if [ "$ok" = 1 ]; then
        read -r K COV NTORD < <("$PYTHON" - "$run" <<'PY'
import json, os, sys
S = json.load(open(os.path.join(sys.argv[1], "summary.json")))
K = S["recommended_K"]
rec = next(r for r in S["frontier"] if r["K"] == K)
print(K, "%.0f%%" % rec["coverage_pct"], rec.get("nt_ordered", rec["nt"]))
PY
)
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" OK "$K" "$COV" "$NTORD" "$run" >> "$SUMMARY"
        echo "  [OK] $name: K=$K, coverage $COV, $NTORD nt ordered"
    else
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" FAILED - - - - >> "$SUMMARY"
        echo "  [FAILED] $name -- see $log"
    fi
    echo
done

echo "=================================================================="
echo "SUMMARY  ($SUMMARY)"
echo "=================================================================="
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
