#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 06:00:00
# 6h walltime (not 36h): each Group has only 2-3 cells and finishes in <1h. A short
# walltime lets the scheduler backfill these into gaps -> starts much sooner on a full partition.
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-64
#SBATCH --job-name=sc2_yap_mc
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2/codes/logs/07.yap_mc/%x.%a.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2/codes/logs/07.yap_mc/%x.%a.err

# Run yap mc mapping, one SLURM array task per snakemake Group command.
# GENOME (hg38 | mm10) is passed via sbatch --export=ALL,GENOME=<genome>.
# snmCseq2 is mixed-species, so there are two independent yap runs:
#   yap_mapping_hg38 (153 cells) and yap_mapping_mm10 (96 cells).
# Mirrors snmCAT/codes/02.run_snakemake.sh (incl. the --rerun-incomplete --unlock fix).

source /home/jmj7858/.bashrc
conda activate mapping
export PATH="$CONDA_PREFIX/bin:$PATH"

: "${GENOME:?Must pass GENOME=hg38 or GENOME=mm10 via sbatch --export}"

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2"
RUN_DIR="${BASE}/yap_mapping_${GENOME}"
CMD_FILE="${RUN_DIR}/snakemake/snakemake_cmd.txt"

mkdir -p "${BASE}/codes/logs/07.yap_mc" "${RUN_DIR}/.slurm_tmp"
cd "${RUN_DIR}/snakemake"

if [[ ! -f "${CMD_FILE}" ]]; then
  echo "ERROR: CMD_FILE not found: ${CMD_FILE}" >&2
  exit 1
fi

N=$(grep -cve '^\s*$' "${CMD_FILE}" || true)
if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${N}" ]]; then
  echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} exceeds number of commands ${N}" >&2
  exit 1
fi

CMD=$(grep -ve '^\s*$' "${CMD_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "[$(date)] GENOME=${GENOME} Task ${SLURM_ARRAY_TASK_ID}/${N}"
echo "CMD: ${CMD}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export TMPDIR="${RUN_DIR}/.slurm_tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"

# Release any stale lock from a cancelled run, then map.
eval "${CMD} --rerun-incomplete --unlock" || true
eval "${CMD} --rerun-incomplete"
echo "[$(date)] DONE"
