#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-15
#SBATCH --job-name=trim.%a
#SBATCH --output=logs/trim/trim.%a.out
#SBATCH --error=logs/trim/trim.%a.err

source ~/.bashrc
conda activate scnomehic

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano

prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' acc_list.txt)

input=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano/fastq
output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano/trimmed_fastq

mkdir -p trimmed_fastq

fastp \
  -i ${input}/${prefix}_1.fq.gz -I ${input}/${prefix}_2.fq.gz \
  -o ${output}/${prefix}_1.t6.fq.gz -O ${output}/${prefix}_2.t6.fq.gz \
  --trim_front1 6 --trim_front2 6 \
  --detect_adapter_for_pe \
  --thread 8 \
  --html ${output}/${prefix}.fastp.html \
  --json ${output}/${prefix}.fastp.json
