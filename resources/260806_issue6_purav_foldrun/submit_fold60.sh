#!/usr/bin/env bash
# Launch the 60% centroid hybrid fold (Part II, Step 15). PRINTS the sbatch commands by default;
# nothing is submitted unless --go is given.
#
#   ./submit_fold60.sh                 # dry run: show what would be submitted
#   TESTN=2 ./submit_fold60.sh --go    # plumbing test: shards 1-2 of each tier only
#   ./submit_fold60.sh --go            # full launch (420 tasks)
#
# Input: data/fold60pid_input/{short_full,mid_fast}.fa from
#   python make_fold_input_hf.py --level 60 --outdir data/fold60pid_input
# Re-running is safe: finished structures are skipped, so resubmitting a tier resumes it.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"; mkdir -p logs

IN=data/fold60pid_input
OUT=data/fold60pid                    # -> $OUT/{short_full,mid_fast}/{structures,metrics}/<NNN>/
MAXCONCURRENT="${MAXCONCURRENT:-128}" # aip-qmorris has no GPU cap; b1 has 168 nodes x 4 L40S
GO=0; [[ "${1:-}" == "--go" ]] && GO=1

declare -A NSHARDS=([short_full]=20 [mid_fast]=400)
for TIER in short_full mid_fast; do
  N=${NSHARDS[$TIER]}
  ARRAY="1-${N}%${MAXCONCURRENT}"; [[ -n "${TESTN:-}" ]] && ARRAY="1-${TESTN}"
  CMD=(sbatch --job-name="fold60-$TIER" --array="$ARRAY"
       --export=ALL,TIER=$TIER,INPUT=$HERE/$IN/$TIER.fa,OUTDIR=$HERE/$OUT/$TIER,NSHARDS=$N,SCRIPTDIR=$HERE
       "$HERE/fold_hybrid.sbatch")
  echo "${CMD[*]}"
  [[ $GO == 1 ]] && "${CMD[@]}"
done
[[ $GO == 1 ]] || echo "(dry run; add --go to submit)"
