#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-51
#SBATCH --job-name=qc
#SBATCH --output=logs/01.qc_trim/qc_trim.%a.txt
#SBATCH --error=logs/01.qc_trim/qc_trim.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1198/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood
conda activate scnomehic

infolder=01.fastq
outfolder=02.fastqc_out
outfolder_fq=03.trimmed_fastq

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}
fastqc --outdir ${outfolder} -t 8 ${infolder}/${prefix}_1.fastq ${infolder}/${prefix}_2.fastq
echo "First fastqc Done"

trim_galore --paired_end -clip_r1 9 -clip_r2 9 -j 8 -o ${outfolder_fq} --gzip --fastqc --fastqc_args "--outdir ${outfolder} -t 8" ${infolder}/${prefix}_1.fastq ${infolder}/${prefix}_2.fastq
echo "Trimming Done"
