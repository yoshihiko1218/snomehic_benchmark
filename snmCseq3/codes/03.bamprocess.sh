#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-100
#SBATCH --job-name=bamprocess
#SBATCH --output=logs/03.bamprocess/bamprocess.%a.out
#SBATCH --error=logs/03.bamprocess/bamprocess.%a.err


source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3
module load java/jdk-17.0.2+8
conda activate scnomehic
# -----------------------
# Paths
# -----------------------
ACC_LIST="acc_list.txt"

vcf=/gpfs/projects/b1198/epifluidlab/shared/data/dbsnp/mm10.dbsnp.vcf
reference=/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10
inputfolder=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3/04.bhmem_bam

prefix="$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' "${ACC_LIST}")"

samtools sort --threads 16 -T ${inputfolder}/${prefix}.tmp -n ${inputfolder}/${prefix}.bhmem.bam | samtools fixmate -m --threads 16 - - | samtools sort --threads 16 -T ${inputfolder}/${prefix}.cor - | samtools markdup -T ${inputfolder}/${prefix}.mdups --threads 16 - - | samtools calmd --threads 16 -b - ${reference}.fa 2>/dev/null > ${inputfolder}/${prefix}.calmd.bam 
samtools index -@ 16 ${inputfolder}/${prefix}.calmd.bam

samtools sort -@ 16 -n -o ${inputfolder}/${prefix}_sorted_by_name.calmd.bam ${inputfolder}/${prefix}.calmd.bam
python /projects/b1198/epifluidlab/yoshii/software/bisulfitehic/src/python/mh_reads_summary.v2.py --in_cram ${inputfolder}/${prefix}_sorted_by_name.calmd.bam --out_summary ${inputfolder}/${prefix}.summary.txt
rm ${inputfolder}/${prefix}_sorted_by_name.calmd.bam

perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 16 --mem 128 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools ${inputfolder}/${prefix}.calmd.bam
perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 16 --mem 128 /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${inputfolder}/${prefix}.calmd.bam ${reference}.fa ${vcf}
