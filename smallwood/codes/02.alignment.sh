#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-51
#SBATCH --job-name=alignment
#SBATCH --output=logs/02.alignment/alignment.%a.txt
#SBATCH --error=logs/02.alignment/alignment.%a.err

# Smallwood 2014 scBS-seq alignment (full workflow, per paper Methods).
#   Step 1  Human-contamination depletion: map trimmed reads to human (hg38) in
#           PAIRED-END mode with Bismark (--bowtie2 --non_directional --unmapped).
#           Reads that do NOT map to human are written to *_unmapped_reads_{1,2}.
#   Step 2  Concatenate the human-unmapped R1+R2 into a single single-end FASTQ
#           (scBS-seq is non-directional; mouse alignment is done single-end).
#   Step 3  Map the human-unmapped reads to mouse (mm10) in SINGLE-END mode
#           (--bowtie2 --non_directional).
# NOTE: paper used GRCh37/NCBIM37; we use hg38/mm10 (mm10 matches the GEO
#       processed *cov.txt build GRCm38). Bismark v0.24.2.

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood

conda activate scnomehic
set -o pipefail

input=03.trimmed_fastq
out_hg38=04.align_hg38
out_mm10=05.align_mm10
ref_hg38=/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/
ref_mm10=/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/

mkdir -p ${out_hg38} ${out_mm10}

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo "[`date`] prefix=${prefix}"

unmap1=${out_hg38}/${prefix}_1_val_1.fq.gz_unmapped_reads_1.fq.gz
unmap2=${out_hg38}/${prefix}_2_val_2.fq.gz_unmapped_reads_2.fq.gz
se_input=${out_hg38}/${prefix}.SE.input.fq.gz
mm10_bam=${out_mm10}/${prefix}.SE.input_bismark_bt2.bam

# ---- Step 1: human (hg38) paired-end alignment, keep unmapped ----
if [ -s "${unmap1}" ] && [ -s "${unmap2}" ]; then
    echo "  [skip] human-unmapped reads already present"
else
    echo "  [run ] Step 1: hg38 PE alignment (--non_directional --unmapped)"
    bismark --bowtie2 --non_directional --unmapped \
        --genome ${ref_hg38} \
        -1 ${input}/${prefix}_1_val_1.fq.gz \
        -2 ${input}/${prefix}_2_val_2.fq.gz \
        -o ${out_hg38}/
fi

# ---- Step 2: concatenate human-unmapped R1+R2 into single-end input ----
if [ -s "${se_input}" ]; then
    echo "  [skip] SE input already present: ${se_input}"
else
    echo "  [run ] Step 2: concatenate unmapped R1+R2 -> ${se_input}"
    zcat ${unmap1} ${unmap2} | gzip > ${se_input}
fi

# ---- Step 3: mouse (mm10) single-end alignment ----
if [ -s "${mm10_bam}" ]; then
    echo "  [skip] mm10 SE BAM already present: ${mm10_bam}"
else
    echo "  [run ] Step 3: mm10 SE alignment (--non_directional)"
    bismark --bowtie2 --non_directional \
        --genome ${ref_mm10} \
        ${se_input} \
        -o ${out_mm10}/
fi

echo "[`date`] DONE prefix=${prefix}"
