#!/usr/bin/env bash
# Optional: archive folded structures to S3 AFTER a batch completes (not during folding).
# Resumable — aws s3 sync only uploads new/changed files, so re-running is cheap and safe.
# Pair with a local cleanup if you're offloading in waves to stay under the disk budget.
#
#   ./sync_to_s3.sh            # dry-run first (prints what would upload)
#   ./sync_to_s3.sh --go       # actually upload
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${OUTDIR:-$HERE/data/fold90pid}"
S3="${S3:-s3://petadex-protein-structures/esmfold2-centroids/}"

FLAGS="--dryrun"
[[ "${1:-}" == "--go" ]] && FLAGS=""

echo "sync $OUTDIR  ->  $S3   (${FLAGS:-live})"
aws s3 sync "$OUTDIR/structures" "${S3}structures/" $FLAGS
aws s3 sync "$OUTDIR/metrics"    "${S3}metrics/"    $FLAGS
