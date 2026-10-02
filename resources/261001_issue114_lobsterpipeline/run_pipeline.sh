#!/usr/bin/env bash
#
# run_pipeline.sh -- LOBSTER sequence-purchasing pipeline, end to end, on a VM.
#
#   preprocess (HMM route) -> design -> validate -> order -> oPools xlsx + figures
#
# INPUT  a CSV whose rows are:   cluster_id,hmm_s3_path   (with a header row)
#        one 90pid cluster per row, paired with the S3 path of ITS HMM.
#
#   ./run_pipeline.sh clusters.csv [out_dir]
#
# WHY ONE CORPUS PASS FOR ALL CLUSTERS
#   fetch_clusters.py streams the ~117 GB ORF corpus once; that pass costs the
#   same for 1 cluster as for 50, so EVERY cluster id is fetched together up front
#   (phase A), then the per-cluster work (phase C) runs on the already-split FASTAs.
#
# REQUIREMENTS (see setup_vm.sh): aws cli, seqkit, hmmer (hmmsearch), python3 with
#   psycopg2 + openpyxl + matplotlib.  Must run in us-east-1 (RDS + S3 in region).
#
# RESUMABLE: a cluster that finished writes clusters/<id>/.done and is skipped on
#   re-run.  A cluster that fails is logged and the batch continues.
#
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
SCR="$HERE/scripts"
PY=${PY:-python3}

CSV=${1:?usage: run_pipeline.sh clusters.csv [out_dir]}
OUT=${2:-$HERE/runs/$(date +%Y%m%d_%H%M%S)}

# ---- design parameters (same for every cluster; override via env) ------------
MAX_LIBRARY=${MAX_LIBRARY:-1000}
MAX_NT=${MAX_NT:-50000}
MAX_OLIGO_NT=${MAX_OLIGO_NT:-330}
K_MAX=${K_MAX:-10}
SEED=${SEED:-0}
MIN_COV=${MIN_COV:-0.80}          # weight_cores truncation cutoff
EVALUE=${EVALUE:-1e-5}            # hmmsearch cutoff when the HMM has no GA line
FIGURES=${FIGURES:-fig2}          # fig2 (fast) | full (re-runs the K sweep) | none
ORF_CORPUS_SIZE=${ORF_CORPUS_SIZE:-117g}

mkdir -p "$OUT"
echo "pipeline out: $OUT"
echo "params: lib<=$MAX_LIBRARY nt<=$MAX_NT oligo<=$MAX_OLIGO_NT K<=$K_MAX seed=$SEED min_cov=$MIN_COV figures=$FIGURES"

# ---- parse the CSV (skip header, strip CR/spaces/blank lines) ----------------
CIDS=(); declare -A HMM_OF
while IFS=, read -r cid hmm _rest; do
    cid=$(echo "$cid" | tr -d ' \r\t'); hmm=$(echo "$hmm" | tr -d ' \r\t')
    [ -z "$cid" ] && continue
    CIDS+=("$cid"); HMM_OF["$cid"]="$hmm"
done < <(tail -n +2 "$CSV")
[ ${#CIDS[@]} -gt 0 ] || { echo "no cluster rows parsed from $CSV (needs a header row)"; exit 1; }
echo "clusters: ${CIDS[*]}"

# =============================================================================
# PHASE A -- one SQL + one 117 GB corpus pass for ALL clusters
# =============================================================================
PRE="$OUT/preprocess"
if [ -f "$PRE/split.done" ]; then
    echo "phase A: already done, skipping corpus pass"
else
    mkdir -p "$PRE"
    echo "phase A: fetch_clusters ids -> extract -> split"
    $PY "$SCR/fetch_clusters.py" ids "${CIDS[@]}" --out-dir "$PRE"
    $PY "$SCR/fetch_clusters.py" extract --out-dir "$PRE" --run
    $PY "$SCR/fetch_clusters.py" split   --out-dir "$PRE"
    touch "$PRE/split.done"
fi

# =============================================================================
# PHASE B -- download each UNIQUE HMM once
# =============================================================================
HMMDIR="$OUT/hmms"; mkdir -p "$HMMDIR"
for hmm in $(printf '%s\n' "${HMM_OF[@]}" | sort -u); do
    dst="$HMMDIR/$(basename "$hmm")"
    [ -f "$dst" ] || { echo "phase B: downloading $hmm"; aws s3 cp "$hmm" "$dst" --no-sign-request; }
done

# =============================================================================
# PHASE C -- per cluster: hmm search -> cores -> weights -> design -> order
# =============================================================================
OK=(); FAILED=()
for cid in "${CIDS[@]}"; do
    cdir="$OUT/clusters/$cid"
    if [ -f "$cdir/.done" ]; then echo "== $cid: already done, skipping"; OK+=("$cid"); continue; fi
    mkdir -p "$cdir"
    echo "== $cid: running (log: $cdir/log.txt)"
    (
        set -e
        foldfa="$PRE/clusters/$cid/cluster.fold.fa"
        [ -f "$foldfa" ] || { echo "no cluster.fold.fa -- cluster $cid absent from corpus?"; exit 2; }
        hmm="$HMMDIR/$(basename "${HMM_OF[$cid]}")"

        if grep -qE '^(GA|TC|NC) ' "$hmm"; then CUT="--cut_ga"; else CUT="-E $EVALUE"; fi
        echo "hmm cutoff: $CUT"

        sed '/^>/!s/\*$//' "$foldfa" > "$cdir/clean.fa"
        hmmsearch $CUT --domtblout "$cdir/$cid.domtbl" -A "$cdir/$cid.sto" \
            "$hmm" "$cdir/clean.fa" > "$cdir/$cid.hmmsearch.txt"

        $PY "$SCR/sto_to_cores.py" "$cdir/$cid.sto" "$cdir/$cid.cores.aln.fa" "$cdir/$cid.cores.fa"
        $PY "$SCR/weight_cores.py" "$cdir/$cid.cores.aln.fa" "$cdir/$cid" --min-cov "$MIN_COV"

        $PY "$SCR/cutsearch_design.py" "$cdir/$cid.weighted.aln.fa" \
            --chemistry gg --cut-search dp \
            --max-library "$MAX_LIBRARY" --max-nt "$MAX_NT" --max-oligo-nt "$MAX_OLIGO_NT" \
            --k-max "$K_MAX" --seed "$SEED" --out-dir "$cdir/design"

        run=$(ls -dt "$cdir/design"/*/ | head -1); run=${run%/}
        echo "$run" > "$cdir/run_dir.txt"

        $PY "$SCR/validate_design.py"     "$run"
        $PY "$SCR/annotate_degenerate.py" "$run"
        $PY "$SCR/make_order.py"          "$run"

        case "$FIGURES" in
            fig2) $PY "$SCR/make_control_figures.py" "$run" --cluster "$cid" --fig2-only --out "$cdir/figures" ;;
            full) $PY "$SCR/make_control_figures.py" "$run" --cluster "$cid"             --out "$cdir/figures" ;;
            none) : ;;
        esac

        touch "$cdir/.done"
    ) > "$cdir/log.txt" 2>&1
    if [ $? -eq 0 ]; then echo "   OK   $cid"; OK+=("$cid"); else echo "   FAIL $cid (see $cdir/log.txt)"; FAILED+=("$cid"); fi
done

# =============================================================================
# PHASE D -- combined oPools order sheet + batch summary (successful clusters)
# =============================================================================
PAIRS=()
for cid in "${OK[@]}"; do
    rd=$(cat "$OUT/clusters/$cid/run_dir.txt" 2>/dev/null) || continue
    [ -f "$rd/order/constructs_simple.csv" ] && PAIRS+=("cluster${cid}=$rd")
done
if [ ${#PAIRS[@]} -gt 0 ]; then
    $PY "$SCR/opools_from_run.py" --out "$OUT/opools_order.xlsx" "${PAIRS[@]}"
fi

# batch summary.tsv from each cluster's summary.json
$PY - "$OUT" "${OK[@]}" <<'PYEOF'
import json, os, sys, glob
out = sys.argv[1]; cids = sys.argv[2:]
rows = [("cluster","rec_K","cores_cov","seqs_cov","library","oligos","nt_ordered")]
for cid in cids:
    rd = open(os.path.join(out,"clusters",cid,"run_dir.txt")).read().strip()
    j = json.load(open(os.path.join(rd,"summary.json")))
    fr = [r for r in j["frontier"] if r.get("K")==j.get("recommended_K")]
    r = fr[0] if fr else {}
    rows.append((cid, j.get("recommended_K"),
                 f"{r.get('n_cores_encoded','?')}", f"{r.get('encoded_weight','?')}",
                 r.get("library","?"), r.get("oligos","?"),
                 r.get("nt_ordered", r.get("nt","?"))))
with open(os.path.join(out,"summary.tsv"),"w") as fh:
    for row in rows: fh.write("\t".join(str(x) for x in row)+"\n")
print("wrote", os.path.join(out,"summary.tsv"))
PYEOF

echo ""
echo "DONE.  ok=${#OK[@]}  failed=${#FAILED[@]}"
[ ${#FAILED[@]} -gt 0 ] && echo "failed clusters: ${FAILED[*]}"
echo "outputs under: $OUT"
