#!/bin/bash
# submit_all.sh - build the manifest, submit the array capped at 20 concurrent.
#
#   ./submit_all.sh                      # everything with a single read pair
#   ./submit_all.sh --skip-done          # only pairs without results
#   ./submit_all.sh --allow-multi-run    # include multi-run pairs too
#   MAXJOBS=10 ./submit_all.sh           # different concurrency cap
#
set -euo pipefail

BASE="${BASE:-/mnt/home/macmaneslab/macmanes/compare}"
MAXJOBS="${MAXJOBS:-20}"
MANIFEST="$BASE/manifest.tsv"

cd "$BASE"
mkdir -p logs results

python build_manifest.py --base "$BASE" --out "$MANIFEST" "$@"

N=$(wc -l < "$MANIFEST")
if [[ "$N" -eq 0 ]]; then
  echo "nothing to submit"
  exit 0
fi

echo
echo "submitting $N jobs, max $MAXJOBS at a time"
sbatch --array="1-${N}%${MAXJOBS}" \
       --export=ALL,BASE="$BASE",MANIFEST="$MANIFEST" \
       run_transrate.sbatch

echo
echo "watch:    squeue -u \$USER"
echo "failures: sacct -j <ARRAYID> --format=JobID,JobName,State,Elapsed,MaxRSS"
echo "collect:  python collect_results.py"
