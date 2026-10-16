#!/usr/bin/env bash
# One-shot health/progress report for the 60% fold. Safe to run any time (read-only).
#   bash fold60_status.sh
cd "$(dirname "$0")" || exit 1
source ../../activate_esmfold2.sh >/dev/null 2>&1
echo "=== fold60 status $(date '+%F %H:%M') ==="

# watcher (detached on klogin02; ps only sees it from that node)
FINDW='ps -u '"$USER"' -o pid=,cmd= | awk '"'"'$2=="bash" && $3 ~ /fold60_watch\.sh$/ {print $1; exit}'"'"''
if [[ "$(hostname -s)" == klogin02 ]]; then W=$(eval "$FINDW")
else W=$(ssh -o BatchMode=yes -o ConnectTimeout=5 klogin02 "$FINDW" 2>/dev/null); fi
if [[ -n "$W" ]]; then echo "watcher: alive (pid $W on klogin02)"; else echo "watcher: NOT RUNNING"; fi
echo "watcher log (last 3):"; tail -3 logs/fold60_watch.log | sed 's/^/  /'

# queue and job outcomes since launch
echo "queue: $(squeue -h -u "$USER" -n fold60-short_full,fold60-mid_fast -o %T | sort | uniq -c | xargs)"
echo "outcomes: $(sacct -X -n -u "$USER" -S 2026-10-15T22:00 --name=fold60-short_full,fold60-mid_fast -o State%12 | sort | uniq -c | xargs)"
NEWFAIL=$(sacct -X -n -u "$USER" -S 2026-10-15T23:53 --name=fold60-short_full,fold60-mid_fast -o JobID%16,State%12,Elapsed,NodeList%10 | awk '$2 ~ /FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL/')
echo "failures/timeouts since 23:53: $( [ -n "$NEWFAIL" ] && echo "$NEWFAIL" | wc -l || echo 0)"
[ -n "$NEWFAIL" ] && echo "$NEWFAIL" | head -10 | sed 's/^/  /'

python - <<'EOF' 2>&1 | grep -v -iE "pynvml|jit|warnings.warn|GPU information"
import csv, glob, os
# progress by predicted fold cost (sweep fits), not by count: long Fast sequences come last
FIT = {"short_full": (0.5021, 0.159), "mid_fast": (8.989e-05, 1.655)}
tot_cost = done_cost = 0.0
for t, (a, b) in FIT.items():
    L = {}
    rid = None
    for line in open(f"data/fold60pid_input/{t}.fa"):
        if line.startswith(">"): rid = line[1:].strip()
        else: L[rid] = len(line.strip())
    done, fails, tmp = set(), {}, 0
    for p in glob.glob(f"data/fold60pid/{t}/progress*.csv"):
        for r in csv.DictReader(open(p)):
            if r["status"] == "ok": done.add(r["id"])
            else: fails[r["status"]] = fails.get(r["status"], 0) + 1
    tmp = len(glob.glob(f"data/fold60pid/{t}/structures/*/*.tmp"))
    c_all = sum(a * l ** b for l in L.values()) / 3600
    c_done = sum(a * L[i] ** b for i in done if i in L) / 3600
    tot_cost += c_all; done_cost += c_done
    print(f"{t:10s}: {len(done):>9,} / {len(L):,} folded ({len(done)/len(L):5.1%}), "
          f"{c_done:6.0f} / {c_all:4.0f} GPU-h of folding ({c_done/c_all:5.1%}); per-seq failures {fails or 0}; in-flight .tmp {tmp}")
print(f"OVERALL   : {done_cost/tot_cost:5.1%} of predicted fold time done ({done_cost:.0f} / {tot_cost:.0f} GPU-h)")
EOF
echo "disk: $(du -sh data/fold60pid 2>/dev/null | cut -f1) | aip-qmorris /project: $(diskusage_report 2>/dev/null | awk '/project aip-qmorris/ {print $5" space, "$6$7" files"; exit}')"
