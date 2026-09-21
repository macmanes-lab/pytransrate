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

# The assembly always gets rewritten into scratch, because TransRate takes the
# contig name from the text before the first '|' or space. ENA writes headers as
#     >ENA|GABB01000001|GABB01000001.1 Species mRNA ...
# so every contig would be named "ENA" and TransRate aborts with
# "Non unique fasta identifier found". NCBI headers (>GABB01000001.1 TSA: ...)
# are already unique; normalising them is harmless.
ASM_LOCAL="$SCRATCH/${PAIR}.fasta"
echo "staging + normalising headers -> $ASM_LOCAL"

if [[ "$ASSEMBLY" == *.gz ]]; then
  READER="gunzip -c"
else
  READER="cat"
fi

$READER "$ASSEMBLY" | awk '
  /^>/ {
    h = substr($0, 2)
    split(h, a, " ")          # drop the description
    id = a[1]
    n = split(id, b, "|")
    if (n >= 3 && b[1] == "ENA") {
      id = b[2]               # ENA|ACCESSION|ACCESSION.V  ->  ACCESSION
    } else if (n > 1) {
      gsub(/\|/, "_", id)     # anything else: keep it all, just make it safe
    }
    print ">" id
    next
  }
  { print }
' > "$ASM_LOCAL"

# refuse to run on a non-unique set rather than let TransRate fail 20 minutes in
total=$(grep -c '^>' "$ASM_LOCAL" || true)
uniq_n=$(grep '^>' "$ASM_LOCAL" | sort -u | wc -l)
echo "contigs: $total, unique names: $uniq_n"
if [[ "$total" -eq 0 ]]; then
  echo "ERROR: no FASTA records in $ASSEMBLY" >&2
  exit 1
fi
if [[ "$total" -ne "$uniq_n" ]]; then
  echo "ERROR: $((total - uniq_n)) duplicate contig names after normalisation." >&2
  echo "       first few:" >&2
  grep '^>' "$ASM_LOCAL" | sort | uniq -d | head -5 >&2
  exit 1
fi
echo "first header: $(grep -m1 '^>' "$ASM_LOCAL")"
echo "scratch used: $(du -sh "$SCRATCH" 2>/dev/null | cut -f1 || echo 0)"

# ---- run
module purge
module load anaconda/colsa
conda activate transrate2026

cd "$RESULTS"

set +e
pytransrate \
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
