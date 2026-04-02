#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-96
#SBATCH --job-name=alignment
#SBATCH --output=logs/02.alignment/alignment.%a.txt
#SBATCH --error=logs/02.alignment/alignment.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/methylhic_new

conda activate scnomehic

module load java/jdk-17.0.2+8

picard=/home/jmj7858/software/picard/picard.jar
reference=/gpfs/projects/b1198/epifluidlab/shared/data/genomes/mm10/mm10
inputfolder=03.trimmed_fastq
outputfolder=04.alignment

mkdir -p ${outputfolder}

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

java -Xmx15G -Djava.library.path=/home/jmj7858/epifluidlab/software/bisulfitehic/jbwa/jbwa-1.0.0/src/main/native -cp "/home/jmj7858/epifluidlab/software/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/home/jmj7858/epifluidlab/software/bisulfitehic/jbwa/jbwa-1.0.0/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem ${reference}.fa -nonDirectional -pbat ${outputfolder}/${prefix}.bam ${inputfolder}/${prefix}_1_val_1.fq.gz ${inputfolder}/${prefix}_2_val_2.fq.gz -t 16 -rgId ${prefix} -rgSm methylHiC -buffer 1000000 -enzymeList /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/dpnII.span_region.bedgraph -outputMateDiffChr

samtools sort --threads 16 -T ${outputfolder}/${prefix}.tmp -n ${outputfolder}/${prefix}.bam | samtools fixmate -m --threads 16 - - | samtools sort --threads 16 -T ${outputfolder}/${prefix}.cor - | samtools markdup -T ${outputfolder}/${prefix}.mdups --threads 16 - - | samtools calmd --threads 16 -b - ${reference}.fa 2>/dev/null > ${outputfolder}/${prefix}.calmd.bam 
samtools index -@ 16 ${outputfolder}/${prefix}.calmd.bam
