#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-15
#SBATCH --job-name=qc.%a
#SBATCH --output=logs/qc/qc.%a.out
#SBATCH --error=logs/qc/qc.%a.err

source ~/.bashrc
conda activate scnomehic

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano

mkdir -p logs/qc
mkdir -p qc_stats

prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' acc_list.txt)

python schic_qc_per_cell.py \
    --cell_id ${prefix} \
    --fastp_dir trimmed_fastq \
    --bam_dir alignment \
    --output_dir qc_stats \
    --mapq 30

