#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-15
#SBATCH --job-name=summary.%a
#SBATCH --output=logs/summary/summary.%a.out
#SBATCH --error=logs/summary/summary.%a.err

source ~/.bashrc
conda activate scnomehic

cd /gpfs/projects/b1198/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano

prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' acc_list.txt)

inputfolder=alignment

samtools sort -@ 16 -n -o ${inputfolder}/${prefix}_sorted_by_name.calmd.bam ${inputfolder}/${prefix}.calmd.bam
python /home/jmj7858/epifluidlab/software/bisulfitehic/src/python/mh_reads_summary.v2.py --in_cram ${inputfolder}/${prefix}_sorted_by_name.calmd.bam --out_summary ${inputfolder}/${prefix}.summary.txt

