#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=scnome_qc_collect
#SBATCH --output=logs/03.qc/qc_all.out
#SBATCH --error=logs/03.qc/qc_all.err

source /home/jmj7858/.bashrc
conda activate scnomehic

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

mkdir -p logs/03.qc qc_stats

echo "=== Per-cell QC collection ($(wc -l < acc_list.txt) cells) ==="
while IFS= read -r prefix; do
    [ -z "$prefix" ] && continue
    echo "  $prefix"
    python codes/scnome_qc_per_cell.py \
        --cell_id "$prefix" \
        --project_dir . \
        --output_dir qc_stats \
        --mapq 30
done < acc_list.txt

echo ""
echo "=== Aggregating into scnome_qc_summary.csv ==="
python codes/collect_scnome_qc.py \
    --qc_dir qc_stats \
    --acc_list acc_list.txt \
    --output scnome_qc_summary.csv

echo ""
echo "Done: scnome_qc_summary.csv"
