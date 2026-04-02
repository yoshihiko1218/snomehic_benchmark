#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-15
#SBATCH --job-name=align.%a
#SBATCH --output=logs/align/align.%a.out
#SBATCH --error=logs/align/align.%a.err

source ~/.bashrc
conda activate bwa

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano
prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' acc_list.txt)

REF=/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa
REF_PREFIX=/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10
input=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano/trimmed_fastq
output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/nagano/alignment

mkdir -p ${output}

RG=$(printf "@RG\tID:%s\tSM:%s\tPL:ILLUMINA" "$prefix" "$prefix")

bowtie2 \
  -x ${REF_PREFIX} \
  -1 ${input}/${prefix}_1.t6.fq.gz \
  -2 ${input}/${prefix}_2.t6.fq.gz \
  --end-to-end \
  --very-sensitive \
  -X 2000 \
  -p 16 \
  --rg-id ${prefix} \
  --rg SM:${prefix} \
  --rg PL:ILLUMINA \
| samtools view -@ 8 -bS - \
| samtools sort -@ 16 -o ${output}/${prefix}.sorted.bam -

samtools index -@ 16 ${output}/${prefix}.sorted.bam

samtools index "${output}/${prefix}.sorted.bam"

samtools sort --threads 16 -T ${output}/${prefix}.tmp -n ${output}/${prefix}.sorted.bam | samtools fixmate -m --threads 16 - - | samtools sort --threads 16 -T ${output}/${prefix}.cor - | samtools markdup -T ${output}/${prefix}.mdups --threads 16 - - | samtools calmd --threads 16 -b - ${REF} 2>/dev/null > ${output}/${prefix}.calmd.bam 
samtools index -@ 16 ${output}/${prefix}.calmd.bam
