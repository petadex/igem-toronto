#!/bin/bash
# Final check for the 60% fold (CPU only): per-ORF completeness, GPU-hours from sacct, full QC,
# disk. Writes data/fold60pid/qc/final_check_<label>.txt. Runs unattended.
#   sbatch --export=ALL,LABEL=l40s fold60_final_check.sbatch
#   sbatch --dependency=afterany:<h100 jobs> --export=ALL,LABEL=all fold60_final_check.sbatch
#SBATCH --job-name=fold60-final-check
#SBATCH --account=aip-qmorris
#SBATCH --cpus-per-task=16
#SBATCH --mem=48gb
#SBATCH --time=04:00:00
#SBATCH --output=logs/fold60-final-check-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source ../../activate_esmfold2.sh >/dev/null 2>&1
OUT=data/fold60pid/qc/final_check_${LABEL:-run}.txt; mkdir -p data/fold60pid/qc
{
echo "=== fold60 final check (${LABEL:-run}) $(date '+%F %H:%M') on $(hostname) ==="
echo; echo "--- per-ORF completeness"
python verify_fold60.py
echo; echo "--- SLURM job outcomes and GPU-hours (all fold60 jobs since launch)"
sacct -X -n -u "$USER" -S 2026-10-15T22:00 --name=fold60-short_full,fold60-mid_fast,fold60-tail \
      -o JobName%20,State%12,ElapsedRaw | awk '{n[$1" "$2]++; h[$1" "$2]+=$3/3600; t+=$3/3600}
      END{for(k in h) printf "  %-34s %5d jobs %8.1f GPU-h\n", k, n[k], h[k]; printf "  TOTAL %41.1f GPU-h\n", t}'
echo; echo "--- full QC (1000 per tier)"
python qc_structures.py --run data/fold60pid --input data/fold60pid_input --per-tier 1000 --workers 16 --seed 99 \
       --out data/fold60pid/qc/final_${LABEL:-run} 2>&1 | grep -v -iE "pynvml|jit|warnings.warn|GPU information"
echo; echo "--- disk"
du -sh data/fold60pid/short_full data/fold60pid/mid_fast data/fold60pid/tail_h100 2>/dev/null
diskusage_report 2>/dev/null | awk '/project aip-qmorris/ {print "  aip-qmorris /project:", $4, "space,", $5, $6, "files"; exit}'
echo "=== done $(date '+%F %H:%M') ==="
} > "$OUT" 2>&1
echo "wrote $OUT"
