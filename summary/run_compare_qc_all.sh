#!/bin/bash
# Run unified QC comparison across all 9 benchmark methods.
# Edit paths to match where MappingSummary.csv.gz / qc_summary.csv files live.

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/summary

source /home/jmj7858/.bashrc
conda activate scnomehic

python compare_qc_all.py \
    --config datasets_all.csv \
    --output all_methods_qc
