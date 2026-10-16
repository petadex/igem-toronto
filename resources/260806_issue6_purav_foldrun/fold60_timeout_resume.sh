#!/usr/bin/env bash
# Resume Fast shards that hit the 3 h limit. Shards 1-5 hold the most short (150-200 aa) sequences
# and need 3.0-3.7 h at the measured 0.74-0.8 s/fold (the power-law cost fit underestimated the
# fixed per-fold overhead at short lengths). Waits for each listed job to END, and only if it ended
# in TIMEOUT resubmits that one shard (fold.py skips finished structures, so it resumes).
#   bash fold60_timeout_resume.sh 6293834_1:1 6293834_2:2 6294471_3:3 6295015_4:4 6295015_5:5
cd "$(dirname "$0")" || exit 1
H=$PWD
declare -A PEND=()
for x in "$@"; do PEND[${x%%:*}]=${x##*:}; done
while ((${#PEND[@]})); do
  for j in "${!PEND[@]}"; do
    if ! squeue -h -j "$j" 2>/dev/null | grep -q .; then
      st=$(sacct -X -n -j "$j" -o State%12 | awk '{print $1}' | tail -1)
      # accounting lags the queue: decide only once sacct shows a terminal state (bug on 02:02:
      # shards 1-2 read as RUNNING right after timing out and were skipped)
      [[ "$st" =~ ^(RUNNING|PENDING|COMPLETING|REQUEUED|CONFIGURING|SUSPENDED)?$ ]] && continue
      k=${PEND[$j]}
      if [[ "$st" == TIMEOUT ]]; then
        new=$(sbatch --parsable --exclude=kn140 --job-name=fold60-mid_fast --array="$k" \
              --export=ALL,TIER=mid_fast,INPUT=$H/data/fold60pid_input/mid_fast.fa,OUTDIR=$H/data/fold60pid/mid_fast,NSHARDS=400,SCRIPTDIR=$H \
              fold_hybrid.sbatch)
        echo "$(date +%H:%M) $j (mid_fast shard $k) TIMEOUT -> resume as $new"
      else
        echo "$(date +%H:%M) $j (mid_fast shard $k) ended $st, no action"
      fi
      unset "PEND[$j]"
    fi
  done
  sleep 60
done
echo "$(date +%H:%M) timeout-resume helper done"
