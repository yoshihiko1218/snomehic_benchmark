#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH -N 1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=scnome_hcg_pc
#SBATCH --output=logs/hcg_genome/scnome_hcg_percell.txt
#SBATCH --error=logs/hcg_genome/scnome_hcg_percell.err

# Per-cell (R1 union R2) destranded HCG for scnome -- correct merge of the two mates.
source /home/jmj7858/.bashrc
conda activate scnomehic
B=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
cd $B/summary
mkdir -p logs/hcg_genome
echo "[$(date)] scnome_hcg_percell.py"
python3 scnome_hcg_percell.py
echo "[$(date)] DONE"
