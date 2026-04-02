#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-41
#SBATCH --job-name=scnome_qc
#SBATCH --output=logs/03.qc/qc.%a.out
#SBATCH --error=logs/03.qc/qc.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

conda activate scnomehic

mkdir -p logs/03.qc
mkdir -p qc_stats

prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' acc_list.txt)
echo "Cell: ${prefix}"

python codes/scnome_qc_per_cell.py \
    --cell_id "${prefix}" \
    --project_dir . \
    --output_dir qc_stats \
    --mapq 30
