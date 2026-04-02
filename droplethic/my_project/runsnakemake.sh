#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=runsnakemake
#SBATCH --output=logs/runsnakemake/runsnakemake.txt
#SBATCH --error=logs/runsnakemake/runsnakemake.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/droplethic/my_project

conda activate rupture

rupture run --configfile config.yaml --cores 16 --snakemake-args --rerun-incomplete --keep-going

