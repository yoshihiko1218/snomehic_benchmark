#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-64
#SBATCH --job-name=snmCAT_mct
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/02.mapping/snakemake.%a.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/02.mapping/snakemake.%a.err

# Run yap mct mapping: one SLURM array task per snakemake Group command.
# Mirrors snmCseq3/codes/01.excute_snakemake.sh.

source /home/jmj7858/.bashrc
conda activate mapping
export PATH="$CONDA_PREFIX/bin:$PATH"

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
CMD_FILE="${BASE}/mapping/snakemake/snakemake_cmd.txt"

mkdir -p "${BASE}/codes/logs/02.mapping" "${BASE}/mapping/.slurm_tmp"
cd "${BASE}/mapping/snakemake"

# ---- sanity checks ----
if [[ ! -f "${CMD_FILE}" ]]; then
  echo "ERROR: CMD_FILE not found: ${CMD_FILE}" >&2
  exit 1
fi

N=$(grep -cve '^\s*$' "${CMD_FILE}" || true)
if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${N}" ]]; then
  echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} exceeds number of commands ${N}" >&2
  exit 1
fi

# Pull the Nth non-empty line (1-indexed)
CMD=$(grep -ve '^\s*$' "${CMD_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}/${N}"
echo "CMD: ${CMD}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export TMPDIR="${BASE}/mapping/.slurm_tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"

# Release any stale lock left by a previously cancelled run, then map.
# --rerun-incomplete redoes outputs that were killed mid-write while reusing valid ones.
eval "${CMD} --unlock" || true
eval "${CMD} --rerun-incomplete"
echo "[$(date)] DONE"
