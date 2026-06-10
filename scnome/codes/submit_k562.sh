#!/bin/bash
# Submit the K562 per-cell pipeline as a dependency chain (11 cells = consecutive
# SRR pairs). Reuses existing per-run bismark BAMs (trim+align already done), so
# the chain is: merge+dedup (per cell) -> NOMe methylation extraction (per cell).
# QC (run_qc_k562) is appended separately once the QC step is built.
set -euo pipefail

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

j1=$(sbatch --parsable codes/03k.merge_dedup_k562.sh)
echo "03k.merge_dedup_k562.sh  -> ${j1}"

j2=$(sbatch --parsable --dependency=afterok:${j1} codes/03k.methy_extract_k562.sh)
echo "03k.methy_extract_k562.sh -> ${j2}  (afterok ${j1})"

echo ""
echo "Submitted K562 chain: ${j1} -> ${j2}"
echo "Monitor: squeue -u jmj7858"
