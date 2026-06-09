#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-12
#SBATCH --job-name=methyext
#SBATCH --output=logs/03.methy_extract/methyext.%a.txt
#SBATCH --error=logs/03.methy_extract/methyext.%a.err

# scNOMe-seq methylation extraction (Pott 2017 protocol), Bismark route.
#
#   Step 1  bismark_methylation_extractor -s --ignore 6 --bedGraph --CX
#           --genome_folder <hg38_bismark>  <rmdup.bam>
#           -> CX coverage of COVERED cytosines (*.bismark.cov.gz)  [per-locus,
#              per-strand, NO SNP filtering -- the YAP/allcools-consistent
#              convention]. NOTE: we use --bedGraph --CX rather than the
#              protocol's --cytosine_report, which would emit a genome-wide
#              report of every cytosine in hg38 (~1.1e9 rows) per cell. The
#              .cov.gz holds the same per-locus calls for covered positions,
#              which is all the NOMe split (step 2) needs.
#
#   Step 2  coverage2cytosine --nome-seq --genome_folder <hg38_bismark>  (input
#           cov must be CX; the --CX flag itself is NOT passed here)
#           -> splits into NOMe CpG (ACG/TCG) and GpC (GCA/GCC/GCT) reports,
#              automatically dropping ambiguous GCG positions. --nome-seq sets the
#              coverage threshold to 1, so only covered positions are reported.
#
# Each read of a pair was aligned independently (mates _1 and _2), so both are
# extracted separately, matching how the libraries were aligned.

source /home/jmj7858/.bashrc
conda activate scnomehic
set -o pipefail

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome
mkdir -p logs/03.methy_extract

BIS=/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/
ALIGN=04.alignment
METHY=05.methy
mkdir -p ${METHY}

# acc_list.txt now holds 34 cells (controls removed): GM12878 = lines 1-12,
# K562 = lines 13-34. Array is scoped to 1-12 = GM12878 (SRR3729642-3729653),
# which were re-trimmed with the both-ends clip. For all cells, set --array=1-34.
prefix=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n' acc_list.txt)
echo "[$(date)] prefix=${prefix}"

for mate in 1 2; do
    bam=${ALIGN}/${prefix}_${mate}.rmdup.bam
    if [ ! -s "${bam}" ]; then
        echo "  WARN: missing ${bam}, skipping"
        continue
    fi

    cov=${METHY}/${prefix}_${mate}.rmdup.bismark.cov.gz

    # ---- Step 1: genome-wide CX coverage (resume: skip if cov already present) ----
    if [ -s "${cov}" ]; then
        echo "  [skip] extractor done: ${cov}"
    else
        echo "  [run ] bismark_methylation_extractor ${bam}"
        bismark_methylation_extractor \
            -s --ignore 6 --comprehensive --multicore 4 \
            --bedGraph --CX \
            --genome_folder "${BIS}" \
            -o "${METHY}" \
            "${bam}"
    fi

    # ---- Step 2: NOMe CpG / GpC split with GCG removed (resume-aware) ----
    cpg=${METHY}/${prefix}_${mate}.NOMe.CpG_report.txt.gz
    gpc=${METHY}/${prefix}_${mate}.NOMe.GpC_report.txt.gz
    if [ -s "${cpg}" ] && [ -s "${gpc}" ]; then
        echo "  [skip] NOMe split done: ${cpg}"
    elif [ -s "${cov}" ]; then
        echo "  [run ] coverage2cytosine --nome-seq ${cov}"
        coverage2cytosine \
            --nome-seq \
            --genome_folder "${BIS}" \
            --dir "${METHY}" \
            --gzip \
            -o "${prefix}_${mate}" \
            "${cov}"
    else
        echo "  WARN: no cov for ${prefix}_${mate}, cannot run coverage2cytosine"
    fi
done

echo "[$(date)] DONE prefix=${prefix}"
