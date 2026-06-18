#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH -N 1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=scnome_loci_pc
#SBATCH --output=logs/hcg_genome/scnome_loci_percell.txt
#SBATCH --error=logs/hcg_genome/scnome_loci_percell.err
source /home/jmj7858/.bashrc
conda activate scnomehic
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/summary
mkdir -p logs/hcg_genome
echo "[$(date)] scnome_loci_percell.py"
python3 scnome_loci_percell.py
echo "[$(date)] DONE"
