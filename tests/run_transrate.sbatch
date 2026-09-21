#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J transrate
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH --exclude=node117,node118
#SBATCH --output=logs/slurm_%A_%a.out
#
# One TransRate run per line of manifest.tsv. Concurrency is capped by the %N
# suffix on the array range, not in here:
#
#     sbatch --array=1-47%20 run_transrate.sbatch
#
# Adjust BASE if the tree lives elsewhere.

set -euo pipefail

BASE="${BASE:-/mnt/home/macmaneslab/macmanes/compare}"
MANIFEST="${MANIFEST:-$BASE/manifest.tsv}"
RESULTS="$BASE/results"
THREADS="${SLURM_CPUS_PER_TASK:-24}"

mkdir -p "$RESULTS" logs

# ---- pick this task's line
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
if [[ -z "$LINE" ]]; then
  echo "no line ${SLURM_ARRAY_TASK_ID} in $MANIFEST" >&2
  exit 1
fi
IFS=$'\t' read -r PAIR ASSEMBLY LEFT RIGHT <<< "$LINE"

# per-pair log, as well as the array log above
exec > >(tee -a "logs/${PAIR}.log") 2>&1

echo "=== $PAIR ==="
echo "node      : $(hostname)"
echo "job       : ${SLURM_ARRAY_JOB_ID:-na}_${SLURM_ARRAY_TASK_ID:-na}"
echo "started   : $(date)"
echo "assembly  : $ASSEMBLY"
echo "read pairs: $(awk -F, '{print NF}' <<< "$LEFT")"

# ---- scratch: TransRate needs the assembly unpacked. Reads stay gzipped.
SCRATCH="${TMPDIR:-/tmp}/transrate.${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH"
cleanup() { rm -rf "$SCRATCH"; }
trap cleanup EXIT

if [[ "$ASSEMBLY" == *.gz ]]; then
  ASM_LOCAL="$SCRATCH/$(basename "${ASSEMBLY%.gz}")"
  echo "unpacking assembly -> $ASM_LOCAL"
  gunzip -c "$ASSEMBLY" > "$ASM_LOCAL"
else
  ASM_LOCAL="$ASSEMBLY"
fi
echo "scratch used: $(du -sh "$SCRATCH" 2>/dev/null | cut -f1 || echo 0)"

# ---- run
module purge
module load anaconda/colsa
conda activate orp

cd "$RESULTS"

set +e
"$HOME"/Oyster_River_Protocol/software/orp-transrate/transrate \
  --assembly "$ASM_LOCAL" \
  --left "$LEFT" \
  --right "$RIGHT" \
  --threads "$THREADS" \
  --output "$RESULTS/$PAIR"
status=$?
set -e
echo "finished  : $(date)  (exit $status)"

if [[ -s "$RESULTS/$PAIR/assemblies.csv" ]]; then
  echo "--- assemblies.csv ---"
  cat "$RESULTS/$PAIR/assemblies.csv"
else
  echo "WARNING: no assemblies.csv written for $PAIR" >&2
  exit 1
fi
