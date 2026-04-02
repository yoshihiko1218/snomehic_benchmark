#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-200
#SBATCH --job-name=rename_fastq_headers
#SBATCH --output=logs/00.rename_fastq/rename_fastq.%a.out
#SBATCH --error=logs/00.rename_fastq/rename_fastq.%a.err

source /home/jmj7858/.bashrc

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3"
IN_DIR="${BASE}/test_fastq"
OUT_DIR="${BASE}/test_fastq_renamed"
SCRIPT="${BASE}/codes/rename_fastq.sh"

mkdir -p "$OUT_DIR" "${BASE}/logs/00.rename_fastq"

cd "$BASE"

mapfile -t FILES < <(find "$IN_DIR" -maxdepth 1 -type f \( -name "*-R1.fq.gz" -o -name "*-R2.fq.gz" \) | sort)

N="${#FILES[@]}"
if [[ "$N" -lt 1 ]]; then
    echo "ERROR: no input files found in $IN_DIR" >&2
    exit 1
fi

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ "$IDX" -lt 0 || "$IDX" -ge "$N" ]]; then
    echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} out of range (1..$N)" >&2
    exit 1
fi

INFILE="${FILES[$IDX]}"
BASENAME="$(basename "$INFILE")"
OUTFILE="${OUT_DIR}/${BASENAME}"

if [[ "$BASENAME" == *-R1.fq.gz ]]; then
    TAG=1
elif [[ "$BASENAME" == *-R2.fq.gz ]]; then
    TAG=2
else
    echo "ERROR: cannot determine read tag from filename: $BASENAME" >&2
    exit 1
fi

echo "[INFO] Task ${SLURM_ARRAY_TASK_ID}/${N}"
echo "[INFO] INFILE : $INFILE"
echo "[INFO] OUTFILE: $OUTFILE"
echo "[INFO] TAG    : $TAG"

lockdir="${OUTFILE}.lock"
if ! mkdir "$lockdir" 2>/dev/null; then
    echo "ERROR: lock exists for $OUTFILE" >&2
    exit 1
fi
trap 'rmdir "$lockdir" >/dev/null 2>&1 || true' EXIT

bash "$SCRIPT" "$INFILE" "$OUTFILE" "$TAG"