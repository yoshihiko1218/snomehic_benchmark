#!/bin/bash
# Create yap-compatible symlinks: SRRxxx_1.fastq.gz -> SRRxxx-R1.fq.gz (and _2 -> -R2)
# yap start-from-cell-fastq needs the pattern "*-R[12].fq.gz"; cell_id = filename prefix before -R1.
set -euo pipefail

FASTQ_DIR="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/fastq"
cd "${FASTQ_DIR}"

n=0
for f1 in *_1.fastq.gz; do
    prefix="${f1%_1.fastq.gz}"
    f2="${prefix}_2.fastq.gz"
    if [[ ! -s "${f2}" ]]; then
        echo "WARN: missing mate for ${prefix}, skipping" >&2
        continue
    fi
    ln -sf "${f1}" "${prefix}-R1.fq.gz"
    ln -sf "${f2}" "${prefix}-R2.fq.gz"
    n=$((n+1))
done
echo "Linked ${n} cells -> *-R[12].fq.gz"
