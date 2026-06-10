#!/bin/bash
# Quick test: run coverage2cytosine --nome-seq on ONE Smallwood cell's CpG-only
# cov to confirm it produces NOMe.CpG.cov.gz (HCG) correctly from a non-CX cov.
source /home/jmj7858/.bashrc
conda activate scnomehic
set -o pipefail

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood
BIS=/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/
cell=SRR1248457
cov=06.methy/${cell}.dedup.bismark.cov.gz
out=06.methy/_hcg_test
mkdir -p "${out}"

echo "[`date`] input cov rows (all-CpG): $(zcat ${cov} | wc -l)"
echo "[`date`] running coverage2cytosine --nome-seq ..."
coverage2cytosine \
    --nome-seq \
    --genome_folder "${BIS}" \
    --dir "${out}" \
    --gzip \
    -o "${cell}" \
    "${cov}"

echo "[`date`] output files:"
ls -la "${out}"
echo "=== HCG (NOMe.CpG) cov rows ==="
zcat "${out}/${cell}.NOMe.CpG.cov.gz" 2>/dev/null | wc -l
echo "=== head of HCG cov ==="
zcat "${out}/${cell}.NOMe.CpG.cov.gz" 2>/dev/null | head -3
echo "=== GpC cov rows (expected ~0 from CpG-only input) ==="
zcat "${out}/${cell}.NOMe.GpC.cov.gz" 2>/dev/null | wc -l
echo "[`date`] DONE"
