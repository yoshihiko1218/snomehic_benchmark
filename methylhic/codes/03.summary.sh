#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-59
#SBATCH --job-name=summary
#SBATCH --output=logs/03.summary/summary.%a.txt
#SBATCH --error=logs/03.summary/summary.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/methylhic

conda activate scnomehic

module load java/jdk-17.0.2+8

picard=/home/jmj7858/software/picard/picard.jar
reference=/gpfs/projects/b1198/epifluidlab/shared/data/genomes/mm10/mm10
vcf=/gpfs/projects/b1198/epifluidlab/shared/data/dbsnp/mm10.dbsnp.vcf
inputfolder=03.trimmed_fastq
outputfolder=04.alignment

mkdir -p ${outputfolder}

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

samtools sort -@ 16 -n -o ${outputfolder}/${prefix}_sorted.calmd.bam ${outputfolder}/${prefix}.calmd.bam
python /gpfs/projects/b1198/epifluidlab/yoshii/software/bisulfitehic/src/python/mh_reads_summary.v2.py --in_cram ${outputfolder}/${prefix}_sorted.calmd.bam --out_summary ${outputfolder}/${prefix}.summary.txt

perl /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 16 --mem 128 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools ${outputfolder}/${prefix}.calmd.bam
perl /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --bistools_path /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 16 --mem 128 /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${outputfolder}/${prefix}.calmd.bam ${reference}.fa ${vcf}
