#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-51
#SBATCH --job-name=bisqc
#SBATCH --output=logs/03.bisqc/bisqc.%a.txt
#SBATCH --error=logs/03.bisqc/bisqc.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood/

conda activate scnomehic

outdir=05.align_mm10
vcf=/gpfs/projects/b1198/epifluidlab/shared/data/dbsnp/mm10.dbsnp.vcf
reference=/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

samtools sort -@ 8 \
  ${outdir}/${prefix}.SE.input_bismark_bt2.bam \
| samtools markdup -@ 8 - \
  ${outdir}/${prefix}.rmdup.bam

RGSM="smallwood"
RGPL="ILLUMINA"
RGLB="lib1"
RGPU="unit1"

samtools addreplacerg \
  -r "@RG\tID:${prefix}\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
  -o "${outdir}/${prefix}.rmdup.RG.bam" \
  "${outdir}/${prefix}.rmdup.bam"

samtools index "${outdir}/${prefix}.rmdup.RG.bam"

perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 8 --mem 64 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools ${outdir}/${prefix}.rmdup.RG.bam

perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 10 --mem 120 /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${outdir}/${prefix}.rmdup.RG.bam ${reference}.fa ${vcf}
