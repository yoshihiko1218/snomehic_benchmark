#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=snmcseq2_qc_collect
#SBATCH --output=logs/04.qc/qc_all.out
#SBATCH --error=logs/04.qc/qc_all.err

source /home/jmj7858/.bashrc
conda activate scnomehic

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2

mkdir -p logs/04.qc qc_stats

echo "=== Per-cell QC collection ($(wc -l < acc_list.txt) cells) ==="
task_id=0
while IFS= read -r prefix; do
    [ -z "$prefix" ] && continue
    task_id=$((task_id + 1))

    # Mirror the genome assignment from 02.alignment.sh
    if [ "$task_id" -le 153 ]; then
        genome="hg38"
    else
        genome="mm10"
    fi

    echo "  [$task_id] $prefix  genome=$genome"
    python codes/snmcseq2_qc_per_cell.py \
        --cell_id "$prefix" \
        --project_dir . \
        --output_dir qc_stats \
        --genome "$genome" \
        --mapq 30
done < acc_list.txt

echo ""
echo "=== Aggregating into snmcseq2_qc_summary.csv ==="
python codes/collect_snmcseq2_qc.py \
    --qc_dir qc_stats \
    --acc_list acc_list.txt \
    --output snmcseq2_qc_summary.csv

echo ""
echo "Done: snmcseq2_qc_summary.csv"
