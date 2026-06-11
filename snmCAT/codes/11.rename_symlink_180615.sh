#!/bin/bash
# yap-pattern symlinks for the 180615 NOMe batch: SRRxxx_1.fastq.gz -> SRRxxx-R1.fq.gz
set -euo pipefail
FASTQ_DIR="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/fastq_180615"
cd "${FASTQ_DIR}"
n=0
for f1 in *_1.fastq.gz; do
    prefix="${f1%_1.fastq.gz}"
    f2="${prefix}_2.fastq.gz"
    [[ -s "${f2}" ]] || { echo "WARN: missing mate ${prefix}" >&2; continue; }
    ln -sf "${f1}" "${prefix}-R1.fq.gz"
    ln -sf "${f2}" "${prefix}-R2.fq.gz"
    n=$((n+1))
done
echo "Linked ${n} cells -> *-R[12].fq.gz in ${FASTQ_DIR}"
