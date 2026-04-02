#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-51
#SBATCH --job-name=alignment
#SBATCH --output=logs/02.alignment/alignment.%a.txt
#SBATCH --error=logs/02.alignment/alignment.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood

conda activate scnomehic

input=03.trimmed_fastq
out_hg38=04.align_hg38
out_mm10=05.align_mm10

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

# zcat ${out_hg38}/${prefix}_1_val_1.fq.gz_unmapped_reads_1.fq.gz \
#     ${out_hg38}/${prefix}_2_val_2.fq.gz_unmapped_reads_2.fq.gz \
# | gzip > ${out_hg38}/${prefix}.SE.input.fq.gz

bismark -bowtie2 -non_directional --genome /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/ ${out_hg38}/${prefix}.SE.input.fq.gz -o ${out_mm10}/ 

# bismark -bowtie2 -non_directional --genome /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/ -1 ${out_hg38}/${prefix}_1_val_1.fq.gz_unmapped_reads_1.fq.gz -2 ${out_hg38}/${prefix}_2_val_2.fq.gz_unmapped_reads_2.fq.gz -o ${out_mm10}/ 
