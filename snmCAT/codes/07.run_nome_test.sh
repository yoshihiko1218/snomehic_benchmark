#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-3
#SBATCH --job-name=nome_test
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/07.nome_test/nt.%a.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/07.nome_test/nt.%a.err

# Re-run 3 cells through yap mct WITH --nome (num_upstr_bases=1 + select-dna-reads --nome)
# to test whether this library has a real NOMe/GCH accessibility signal.
source /home/jmj7858/.bashrc
conda activate mapping
export PATH="$CONDA_PREFIX/bin:$PATH"

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
CMD_FILE="${BASE}/nome_test/mapping/snakemake/snakemake_cmd.txt"
mkdir -p "${BASE}/codes/logs/07.nome_test" "${BASE}/nome_test/.slurm_tmp"

CMD=$(grep -ve '^\s*$' "${CMD_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ${CMD}"
export TMPDIR="${BASE}/nome_test/.slurm_tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"
eval "${CMD} --rerun-incomplete --unlock" || true
eval "${CMD} --rerun-incomplete"
echo "[$(date)] DONE"
