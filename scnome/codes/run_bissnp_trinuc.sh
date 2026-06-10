#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=bistrinuc
#SBATCH --output=logs/05.bissnp_trinuc/trinuc.%a.out
#SBATCH --error=logs/05.bissnp_trinuc/trinuc.%a.err

# BisSNP trinucleotide methylation QC (the ORIGINAL Bis-tools method) on chrM and
# chr21, per cell per mate. Runs BisulfiteGenotyper (-sm BM, EMIT_ALL_CYTOSINES,
# -minPatConv 0.8, -minConv 1, dbSNP-aware, -nonDirectional) over all 16 N-C-N
# contexts, then extracts the per-trinucleotide methylation summary into
# 04.alignment/<cell>_<mate>.rmdup.RG.trinuc_methy.{chrM,chr21}.txt -- exactly the
# files scnome_qc_per_cell.py expects (ACT=conversion control, ACG=CpG, GCT=GpC).
#
# BisSNP is GATK-3.8-based and its walker discovery FAILS on Java >8; must use
# Java 8 (module java/jdk1.8.0_191), NOT the conda env's Java 21.
#
#   GM   : sbatch --array=1-12 --export=ALL,LISTFILE=acc_list.txt             codes/run_bissnp_trinuc.sh
#   K562 : sbatch --array=1-11 --export=ALL,LISTFILE=acc_list_k562_cellids.txt codes/run_bissnp_trinuc.sh

source /home/jmj7858/.bashrc 2>/dev/null
module load java/jdk1.8.0_191
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome || { echo "ERROR: cd failed"; exit 1; }
mkdir -p logs/05.bissnp_trinuc

export BISTOOLS=/projects/b1198/epifluidlab/yoshii/software/Bis-tools
TRINUC=${BISTOOLS}/Bis-QC/after_reads_mapping/bissnp_trinuc_sample.pl
JAR=${BISTOOLS}/Bis-SNP/Bis-SNP.latest.jar
GENOME=/projects/b1198/epifluidlab/yoshii/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
DBSNP=/projects/b1198/epifluidlab/yoshii/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf
ALIGN=04.alignment
LISTFILE=${LISTFILE:-acc_list.txt}

prefix=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n' ${LISTFILE})
echo "[$(date)] cell=${prefix} (list=${LISTFILE})  java=$(which java)"
if [ -z "${prefix}" ]; then echo "ERROR: empty prefix for task ${SLURM_ARRAY_TASK_ID}"; exit 1; fi

for m in 1 2; do
  bam=${ALIGN}/${prefix}_${m}.rmdup.RG.bam
  if [ ! -s "${bam}" ]; then echo "ERROR: missing ${bam}"; exit 1; fi
  for chrom in chrM chr21; do
    out=${ALIGN}/${prefix}_${m}.rmdup.RG.trinuc_methy.${chrom}.txt
    if [ -s "${out}" ]; then echo "[skip] ${out}"; continue; fi
    echo "[run ] BisSNP trinuc ${prefix} mate ${m} ${chrom}"
    perl "${TRINUC}" --bissnp "${JAR}" --genome "${GENOME}" --dbsnp "${DBSNP}" \
        --nt 4 --mem 22 --interval "${chrom}" --nonDirectional "${out}" "${bam}" \
        || { echo "ERROR: BisSNP trinuc failed ${prefix} m${m} ${chrom}"; exit 1; }
    # tidy GATK index left behind for the (already-removed) intermediate vcf
    rm -f "${ALIGN}/${prefix}_${m}.rmdup.RG.cytosine.raw.vcf.idx"
  done
done
echo "[$(date)] DONE ${prefix}"
