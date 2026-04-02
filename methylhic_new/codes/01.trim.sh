#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-96
#SBATCH --job-name=qc
#SBATCH --output=logs/01.qc_trim/qc_trim.%a.txt
#SBATCH --error=logs/01.qc_trim/qc_trim.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/methylhic_new

conda activate scnomehic

infolder=fastq
outfolder=02.fastqc_out
outfolder_fq=03.trimmed_fastq

mkdir -p ${infolder} ${outfolder} ${outfolder_fq}

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}
fastqc --outdir ${outfolder} -t 8 ${infolder}/${prefix}_1.fq.gz ${infolder}/${prefix}_2.fq.gz
echo "First fastqc Done"

trim_galore --paired -j 8 -o ${outfolder_fq} --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --gzip --fastqc --fastqc_args "--outdir ${outfolder} -t 8" ${infolder}/${prefix}_1.fq.gz ${infolder}/${prefix}_2.fq.gz
echo "Trimming Done"
