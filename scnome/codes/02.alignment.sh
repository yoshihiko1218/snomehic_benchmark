#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-41
#SBATCH --job-name=bisqc
#SBATCH --output=logs/02.bisqc/bisqc.%a.txt
#SBATCH --error=logs/02.bisqc/bisqc.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

conda activate scnomehic

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

input=03.trimmed_fastq
outdir=04.alignment
methydir=05.methy
vcf=/projects/b1198/epifluidlab/yoshii/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf
reference=/projects/b1198/epifluidlab/yoshii/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set

# bismark --fastq --prefix ${prefix} --output_dir ${outdir} --non_directional --phred33-quals --score_min L,0,-0.2 --bowtie2 --genome /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/ ${input}/${prefix}_1_trimmed.fq.gz
# bismark --fastq --prefix ${prefix} --output_dir ${outdir} --non_directional --phred33-quals --score_min L,0,-0.2 --bowtie2 --genome /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/ ${input}/${prefix}_2_trimmed.fq.gz

# samtools sort -@ 8 \
#   ${outdir}/${prefix}.${prefix}_1_trimmed_bismark_bt2.bam \
# | samtools markdup -@ 8 - \
#   ${outdir}/${prefix}_1.rmdup.bam

# samtools sort -@ 8 \
#   ${outdir}/${prefix}.${prefix}_2_trimmed_bismark_bt2.bam \
# | samtools markdup -@ 8 - \
#   ${outdir}/${prefix}_2.rmdup.bam


#bismark_methylation_extractor -s --ignore 6 --output ${methydir} --cytosine_report --CX --genome_folder /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/ ${outdir}/${prefix}_1.rmdup.bam
#bismark_methylation_extractor -s --ignore 6 --output ${methydir} --cytosine_report --CX --genome_folder /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/ ${outdir}/${prefix}_2.rmdup.bam

# samtools sort -@ 16 -n -o ${outdir}/${prefix}_1.sorted.rmdup.bam ${outdir}/${prefix}_1.rmdup.bam
# python /gpfs/projects/b1198/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome/codes/se_bam_summary.py --in_bam ${outdir}/${prefix}_1.sorted.rmdup.bam --out_summary ${outdir}/${prefix}_1.summary.txt

# RGSM="scNOMe"
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

# perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 16 --mem 128 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools ${outdir}/${prefix}_1.rmdup.RG.bam
# perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-QC/Bis-QC.pl --QC_mode 1 --disable_enzyme_eff_check --disable_coverage_check --pattern WCH --nt 16 --mem 128 --genome ${reference}.fa --dbsnp ${vcf} --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools ${outdir}/${prefix}_2.rmdup.RG.bam

#rm ${outdir}/${prefix}_1.rmdup.bam
#rm ${outdir}/${prefix}_2.rmdup.bam

# perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --nomeseq --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 10 --mem 120 /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${outdir}/${prefix}_1.rmdup.RG.bam ${reference}.fa ${vcf}
# perl /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates --nomeseq --bistools_path /projects/b1198/epifluidlab/yoshii/software/Bis-tools --lowCov --mmq 30 --nt 10 --mem 120 /projects/b1198/epifluidlab/yoshii/software/Bis-tools/Bis-SNP/BisSNP-0.90.jar ${outdir}/${prefix}_2.rmdup.RG.bam ${reference}.fa ${vcf}


python /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/bam_summary_universal.py --in_bam ${outdir}/${prefix}_1.rmdup.RG.bam  --out_summary  ${outdir}/${prefix}_1.summary.txt
python /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/bam_summary_universal.py --in_bam ${outdir}/${prefix}_2.rmdup.RG.bam  --out_summary  ${outdir}/${prefix}_2.summary.txt