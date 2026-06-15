#!/bin/bash
set -euo pipefail
FASTQ_DIR="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/fastq_brain"
cd "${FASTQ_DIR}"
n=0
for f1 in *_1.fastq.gz; do
    p="${f1%_1.fastq.gz}"; f2="${p}_2.fastq.gz"
    [[ -s "${f2}" ]] || { echo "WARN: missing mate ${p}" >&2; continue; }
    ln -sf "${f1}" "${p}-R1.fq.gz"; ln -sf "${f2}" "${p}-R2.fq.gz"; n=$((n+1))
done
echo "Linked ${n} cells in ${FASTQ_DIR}"
