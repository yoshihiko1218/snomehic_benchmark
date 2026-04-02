#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-250
#SBATCH --job-name=snmcseq2_qc
#SBATCH --output=logs/04.qc/qc.%a.out
#SBATCH --error=logs/04.qc/qc.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2

conda activate scnomehic

mkdir -p logs/04.qc
mkdir -p qc_stats

prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' acc_list.txt)
echo "Cell: ${prefix}"

# Determine genome from task ID (mirrors 02.alignment.sh and 03.summary.sh)
if [ "${SLURM_ARRAY_TASK_ID}" -le 153 ]; then
    genome="hg38"
else
    genome="mm10"
fi
echo "Genome: ${genome}"

python codes/snmcseq2_qc_per_cell.py \
    --cell_id "${prefix}" \
    --project_dir . \
    --output_dir qc_stats \
    --genome "${genome}" \
    --mapq 30
