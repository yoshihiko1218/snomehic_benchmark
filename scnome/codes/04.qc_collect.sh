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

# GM12878 = 12 single-run cells (lines 1-12 of acc_list.txt; prefix = SRR id).
echo "=== Per-cell QC: GM12878 (12 single-run cells) ==="
head -12 acc_list.txt | while IFS= read -r prefix; do
    [ -z "$prefix" ] && continue
    echo "  GM $prefix"
    python codes/scnome_qc_per_cell.py \
        --cell_id "$prefix" --project_dir . --output_dir qc_stats --mapq 30
done

# K562 = 11 merged cells, each = 2 runs (acc_list_k562_cells.tsv: cell run1 run2).
echo "=== Per-cell QC: K562 (11 merged cells, 2 runs each) ==="
while IFS=$'\t' read -r cell r1 r2; do
    [ -z "$cell" ] && continue
    echo "  K562 $cell ($r1,$r2)"
    python codes/scnome_qc_per_cell.py \
        --cell_id "$cell" --runs "${r1},${r2}" \
        --project_dir . --output_dir qc_stats --mapq 30
done < acc_list_k562_cells.tsv

# Combined cell-id list for aggregation (12 GM SRR ids + 11 K562_NN ids).
head -12 acc_list.txt > acc_list_cells_all.txt
cut -f1 acc_list_k562_cells.tsv >> acc_list_cells_all.txt

echo ""
echo "=== Aggregating into scnome_qc_summary.csv (23 cells) ==="
python codes/collect_scnome_qc.py \
    --qc_dir qc_stats \
    --acc_list acc_list_cells_all.txt \
    --output scnome_qc_summary.csv

echo ""
echo "Done: scnome_qc_summary.csv"
