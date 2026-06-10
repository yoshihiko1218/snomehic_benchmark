#!/bin/bash
# Create yap-compatible cell-level FASTQ symlinks for snmC-seq2, SPLIT BY GENOME.
#
# snmCseq2 is a MIXED-SPECIES dataset (Luo et al 2018 snmC-seq2): 153 hg38 cells +
# 96 mm10 cells. snmCseq2_genome_map.tsv maps each mate (SRRxxx_1 / SRRxxx_2) to its
# genome; _1 and _2 of a cell always share the same genome.
#
# yap takes ONE bismark reference per run, so we make two yap runs. This script links:
#   01.fastq/SRRxxx_{1,2}.fq.gz  ->  fastq_yap_<genome>/SRRxxx-R{1,2}.fq.gz
# and writes codes/cells_<genome>.txt (one cell_id per line).
set -euo pipefail

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2"
SRC="${BASE}/01.fastq"
MAP="${BASE}/snmCseq2_genome_map.tsv"

# Build per-cell genome assignment (cell_id = SRR prefix without _1/_2 mate suffix).
declare -A CELL_GENOME
while IFS=$'\t' read -r mate genome; do
    [[ -z "${mate}" ]] && continue
    cell="${mate%_[12]}"
    CELL_GENOME["${cell}"]="${genome}"
done < "${MAP}"

declare -A N
> "${BASE}/codes/cells_hg38.txt"
> "${BASE}/codes/cells_mm10.txt"

for cell in "${!CELL_GENOME[@]}"; do
    genome="${CELL_GENOME[$cell]}"
    f1="${SRC}/${cell}_1.fq.gz"
    f2="${SRC}/${cell}_2.fq.gz"
    if [[ ! -s "${f1}" || ! -s "${f2}" ]]; then
        echo "WARN: missing fastq for ${cell}, skipping" >&2
        continue
    fi
    dst="${BASE}/fastq_yap_${genome}"
    mkdir -p "${dst}"
    ln -sf "${f1}" "${dst}/${cell}-R1.fq.gz"
    ln -sf "${f2}" "${dst}/${cell}-R2.fq.gz"
    echo "${cell}" >> "${BASE}/codes/cells_${genome}.txt"
    N[$genome]=$(( ${N[$genome]:-0} + 1 ))
done

sort -o "${BASE}/codes/cells_hg38.txt" "${BASE}/codes/cells_hg38.txt"
sort -o "${BASE}/codes/cells_mm10.txt" "${BASE}/codes/cells_mm10.txt"

echo "hg38 cells linked: ${N[hg38]:-0}"
echo "mm10 cells linked: ${N[mm10]:-0}"
