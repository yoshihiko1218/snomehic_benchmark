#!/bin/bash
# Compare QC metrics across four datasets: MethylHiC, snm3C-seq3, methylhic_new, scnomehic
# Run from benchmark/ directory, or cd to summary/ and adjust paths.

set -e
cd "$(dirname "$0")"

# Activate conda env with matplotlib/pandas
eval "$(conda shell.bash hook)"
conda activate mapping

python3 compare_qc.py \
  --dataset1 ../methylhic/alignment/stats/MappingSummary.csv.gz \
  --dataset2 ../snmCseq3/alignment/stats/MappingSummary.csv.gz \
  --label1 MethylHiC \
  --label2 snm3C-seq3 \
  --dataset3 ../methylhic_new/alignment/stats/MappingSummary.csv.gz \
  --label3 methylhic_new \
  --dataset4 ../scnomehic/alignment/stats/MappingSummary.csv.gz \
  --label4 scnomehic \
  --output compare_qc_4datasets_methylhic_snm3c_methylhicnew_scnomehic

echo "Done. Outputs:"
echo "  compare_qc_4datasets_methylhic_snm3c_methylhicnew_scnomehic.combined.csv"
echo "  compare_qc_4datasets_methylhic_snm3c_methylhicnew_scnomehic.summary.csv"
echo "  compare_qc_4datasets_methylhic_snm3c_methylhicnew_scnomehic.plots.pdf"
