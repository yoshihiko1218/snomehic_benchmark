#!/bin/bash
# Submit the GM12878 rerun as a dependency chain (array 1-12 = GM cells).
# Each step starts only if the previous array completes successfully (afterok).
# Run AFTER clearing stale GM outputs (codes/util_clear_gm_stale.sh --yes).
set -euo pipefail

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

j1=$(sbatch --parsable codes/01.trim.sh)
echo "01.trim.sh           -> ${j1}"

j2=$(sbatch --parsable --dependency=afterok:${j1} codes/02.align_dedup.sh)
echo "02.align_dedup.sh    -> ${j2}  (afterok ${j1})"

j3=$(sbatch --parsable --dependency=afterok:${j2} codes/03.methy_extract.sh)
echo "03.methy_extract.sh  -> ${j3}  (afterok ${j2})"

j4=$(sbatch --parsable --dependency=afterok:${j3} codes/04.qc_per_cell.sh)
echo "04.qc_per_cell.sh    -> ${j4}  (afterok ${j3})"

echo ""
echo "Submitted GM rerun chain: ${j1} -> ${j2} -> ${j3} -> ${j4}"
echo "Monitor: squeue -u jmj7858"
