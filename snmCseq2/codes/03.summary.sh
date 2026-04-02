#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-250
#SBATCH --job-name=alignment
#SBATCH --output=logs/03.summary/summary.%a.txt
#SBATCH --error=logs/03.summary/summary.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2

conda activate scnomehic

outdir=05.align
outputfolder=05.align
methydir=06.summary
# ==============================
# Dynamic reference selection
# ==============================
if [ ${SLURM_ARRAY_TASK_ID} -le 153 ]; then
    echo "Using HUMAN reference (hg38)"
    vcf=/projects/b1198/epifluidlab/yoshii/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf
    reference=/projects/b1198/epifluidlab/yoshii/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set
else
    echo "Using MOUSE reference (mm10)"
    vcf=/gpfs/projects/b1198/epifluidlab/shared/data/dbsnp/mm10.dbsnp.vcf
    reference=/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10
fi

echo "Reference: ${reference}"
echo "VCF: ${vcf}"
prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

input=05.align
methydir=06.summary

# samtools sort -@ 8 \
#      ${input}/${prefix}_1.clean_bismark_bt2.bam \
#      | samtools markdup -@ 8 - \
#      ${input}/${prefix}_1.rmdup.bam

# samtools sort -@ 8 \
#      ${input}/${prefix}_2.clean_bismark_bt2.bam \
#      | samtools markdup -@ 8 - \
#      ${input}/${prefix}_2.rmdup.bam

# bismark_methylation_extractor -s --ignore 6 --output ${methydir} --cytosine_report --CX --genome_folder /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/ ${input}/${prefix}_1.rmdup.bam
# bismark_methylation_extractor -s --ignore 6 --output ${methydir} --cytosine_report --CX --genome_folder /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/ ${input}/${prefix}_2.rmdup.bam


# RGSM="snmCseq2"
# RGPL="ILLUMINA"
# RGLB="lib1"
# RGPU="unit1"

# samtools addreplacerg \
#   -r "@RG\tID:${prefix}_1\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
#   -o "${outdir}/${prefix}_1.rmdup.RG.bam" \
#   "${outdir}/${prefix}_1.rmdup.bam"

# samtools index "${outdir}/${prefix}_1.rmdup.RG.bam"

# samtools addreplacerg \
#   -r "@RG\tID:${prefix}_2\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
#   -o "${outdir}/${prefix}_2.rmdup.RG.bam" \
#   "${outdir}/${prefix}_2.rmdup.bam"

# samtools index "${outdir}/${prefix}_2.rmdup.RG.bam"


# perl /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 16 --mem 128 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools ${input}/${prefix}_1.rmdup.RG.bam
# perl /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --bistools_path /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 16 --mem 128 /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${input}/${prefix}_1.rmdup.RG.bam ${reference}.fa ${vcf}

# perl /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 16 --mem 128 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools ${input}/${prefix}_2.rmdup.RG.bam
# perl /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --bistools_path /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 16 --mem 128 /gpfs/projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${input}/${prefix}_2.rmdup.RG.bam ${reference}.fa ${vcf}

# samtools addreplacerg \
#   -r "@RG\tID:${prefix}_2\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
#   -o "${outdir}/${prefix}.rmdup.RG.bam" \
#   "${outdir}/${prefix}_1.clean_bismark_bt2_pe.bam"
# samtools sort -@ 8 -o "${outdir}/${prefix}.rmdup.RG.sorted.bam" "${outdir}/${prefix}.rmdup.RG.bam"
# samtools index -@ 8 "${outdir}/${prefix}.rmdup.RG.sorted.bam"
# perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 8 --mem 64 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools ${outdir}/${prefix}.rmdup.RG.sorted.bam

# samtools sort --threads 16 -T ${outputfolder}/${prefix}.tmp -n ${outputfolder}/${prefix}.rmdup.RG.sorted.bam | samtools fixmate -m --threads 16 - - | samtools sort --threads 16 -T ${outputfolder}/${prefix}.cor - | samtools markdup -T ${outputfolder}/${prefix}.mdups --threads 16 - - | samtools calmd --threads 16 -b - ${reference}.fa 2>/dev/null > ${outputfolder}/${prefix}.calmd.bam 
# samtools sort -@ 16 -n -o ${outputfolder}/${prefix}_sorted.calmd.bam ${outputfolder}/${prefix}.calmd.bam
# python /gpfs/projects/b1198/epifluidlab/yoshii/software/bisulfitehic/src/python/mh_reads_summary.v2.py --in_cram ${outputfolder}/${prefix}_sorted.calmd.bam --out_summary ${outputfolder}/${prefix}.summary.txt


# rm ${outdir}/${prefix}_2.rmdup.bam
# rm ${outdir}/${prefix}_1.rmdup.bam
python codes/se_bam_summary.py --in_bam ${outdir}/${prefix}_1.clean_bismark_bt2.bam --out_summary ${outdir}/${prefix}_1.summary.txt 
python codes/se_bam_summary.py --in_bam ${outdir}/${prefix}_2.clean_bismark_bt2.bam --out_summary ${outdir}/${prefix}_2.summary.txt 