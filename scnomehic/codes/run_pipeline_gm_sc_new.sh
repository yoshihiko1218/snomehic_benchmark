#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# From-scratch rerun of the UPDATED sc_NOMeHiC_pipeline on gm_sc_new (hg38).
#
# Scope (user 2026-05-28): trim -> align -> bamprocess -> bisqc -> bistools ->
# methylation. NO hicluster. --forceall regenerates every rule output for all
# 188 cells (alignment + bistools steps were fixed by the user).
#
# This is the DRIVER job: a lightweight long-walltime allocation that runs
# snakemake with the slurm profile, which in turn submits each rule as its own
# child SLURM job (account b1042, partition genomics, <=48h each) via -j 1000.
# Driver lives on genomicslong (10-day max) so it outlives the multi-day run.
#
# Submit:  sbatch codes/run_pipeline_gm_sc_new.sh
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=scnh_driver
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnomehic/logs/pipeline_rerun/driver.%j.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnomehic/logs/pipeline_rerun/driver.%j.err

set -eo pipefail

source /projects/b1198/epifluidlab/yoshii/software/conda/etc/profile.d/conda.sh
conda activate scnomehic
module load java/jdk-17.0.2+8

PIPELINE=/gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline
cd "${PIPELINE}"

echo "=== Driver start: $(date) on $(hostname) ==="
echo "=== sbatch availability check ==="
which sbatch && sbatch --version

# Per-rule SLURM submission via the slurm profile (-j 1000). --forceall makes
# every rule rerun from scratch, overwriting existing outputs in place.
# The Snakefile sets workdir -> gm_sc_new and start_from: raw from configs/config.yaml.
snakemake \
  -s "${PIPELINE}/Snakefile" \
  --profile "${PIPELINE}/profiles/slurm" \
  --configfile "${PIPELINE}/configs/config.yaml" \
  --forceall \
  -j 1000 \
  -p

echo "=== Driver done: $(date) ==="
