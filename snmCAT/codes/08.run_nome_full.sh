#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-64
#SBATCH --job-name=snmCAT_nome
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/08.nome_full/snakemake.%a.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/08.nome_full/snakemake.%a.err

# Full 100-cell re-run of yap mct WITH --nome (num_upstr_bases=1 + select-dna-reads --nome),
# i.e. the paper-correct snmCAT-seq / NOMe processing. Output in mapping_nome/.
source /home/jmj7858/.bashrc
conda activate mapping
export PATH="$CONDA_PREFIX/bin:$PATH"

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
CMD_FILE="${BASE}/mapping_nome/snakemake/snakemake_cmd.txt"
mkdir -p "${BASE}/codes/logs/08.nome_full" "${BASE}/mapping_nome/.slurm_tmp"

CMD=$(grep -ve '^\s*$' "${CMD_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ${CMD}"
export TMPDIR="${BASE}/mapping_nome/.slurm_tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"

eval "${CMD} --rerun-incomplete --unlock" || true
eval "${CMD} --rerun-incomplete"
echo "[$(date)] DONE"
