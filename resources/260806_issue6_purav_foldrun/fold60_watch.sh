#!/usr/bin/env bash
# Babysit the 60% fold (run on the login node, log to logs/fold60_watch.log): (1) quarantine black-hole nodes (task FAILED in < 120 s) by excluding them
# from all pending fold60 arrays and resubmitting the failed task elsewhere; (2) raise the Fast
# throttle in stages so total concurrency stays ~64. Prints one line per action.
cd /project/6122149/purav/igem-toronto/resources/260806_issue6_purav_foldrun || exit 1
H=$PWD
SHORT="6294470 6295014"; TEST="6293833 6293834"; MAIN=6294471; RETRY=6295015
declare -A BAD=([kn140]=1) HANDLED=() RESUB=()
stage=0
START=$(date +%FT%T)
# every failure that already exists has been dealt with (kn140 batch -> 6295014 / 6295015)
while read -r j; do HANDLED[$j]=1; done < <(sacct -X -n -u "$USER" -S "$(date -d '-12 hours' +%FT%T)" \
    --name=fold60-short_full,fold60-mid_fast -o JobID%16,State%10 | awk '$2=="FAILED" {print $1}')
echo "$(date +%H:%M) watcher start: ${#HANDLED[@]} pre-existing failures marked handled"
excl() { local IFS=,; echo "${!BAD[*]}"; }
while :; do
  ALL=$(squeue -h -u "$USER" -o "%A" -n fold60-short_full,fold60-mid_fast | sort -u | paste -sd,)
  # --- black-hole detection over every fold60 array this user owns ---
  while read -r jid state el node; do
    [[ -z "$jid" || -n "${HANDLED[$jid]:-}" ]] && continue
    secs=$(awk -F: '{print ($1*3600)+($2*60)+$3}' <<<"$el")
    if [[ "$state" == FAILED && "$secs" -lt 120 ]]; then
      HANDLED[$jid]=1
      if [[ -z "${BAD[$node]:-}" ]]; then
        BAD[$node]=1
        echo "$(date +%H:%M) quarantine $node (fast failure $jid)"
      fi
      for a in $(squeue -h -u "$USER" -t PD -o "%F" -n fold60-short_full,fold60-mid_fast | sort -u); do
        scontrol update JobId="$a" ExcNodeList="$(excl)" >/dev/null 2>&1
      done
      arr=${jid%_*}; task=${jid#*_}
      name=$(sacct -X -n -j "$jid" -o JobName%30 | xargs)
      tier=${name#fold60-}; n=20; [[ $tier == mid_fast ]] && n=400
      if [[ -n "${RESUB[$tier:$task]:-}" ]]; then echo "$(date +%H:%M) $tier shard $task failed again ($jid), NOT resubmitting"; continue; fi
      RESUB[$tier:$task]=1
      new=$(sbatch --parsable --exclude="$(excl)" --job-name="$name" --array="$task" \
            --export=ALL,TIER=$tier,INPUT=$H/data/fold60pid_input/$tier.fa,OUTDIR=$H/data/fold60pid/$tier,NSHARDS=$n,SCRIPTDIR=$H \
            fold_hybrid.sbatch)
      echo "$(date +%H:%M) resubmitted $jid as $new (exclude $(excl))"
    fi
  done < <(sacct -X -n -u "$USER" -S "$START" --name=fold60-short_full,fold60-mid_fast \
           -o JobID%16,State%10,Elapsed,NodeList%10 2>/dev/null | awk '$2=="FAILED"')
  # --- staged throttle ---
  if [[ $stage == 0 ]] && ! squeue -h -j "${SHORT// /,},${TEST// /,}" 2>/dev/null | grep -q .; then
    scontrol update JobId=$MAIN ArrayTaskThrottle=48 >/dev/null 2>&1
    echo "$(date +%H:%M) short_full + test shards done -> $MAIN throttle 48 (+ $RETRY at 16)"; stage=1
  fi
  if [[ $stage == 1 ]] && ! squeue -h -j $RETRY 2>/dev/null | grep -q .; then
    scontrol update JobId=$MAIN ArrayTaskThrottle=64 >/dev/null 2>&1
    echo "$(date +%H:%M) retry array done -> $MAIN throttle 64"; stage=2
  fi
  # --- done ---
  if ! squeue -h -u "$USER" -n fold60-short_full,fold60-mid_fast 2>/dev/null | grep -q .; then
    echo "$(date +%H:%M) all fold60 arrays finished"; exit 0
  fi
  sleep 120
done
