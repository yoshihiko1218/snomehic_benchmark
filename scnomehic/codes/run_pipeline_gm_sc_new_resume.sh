#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# RESUME run of sc_NOMeHiC_pipeline on gm_sc_new (hg38).
#
# The original from-scratch --forceall run (driver 9439299) completed 182/188
# cells. The 6 LARGEST cells repeatedly TIMED OUT in demultiplex_fastqc_trim:
# the rule only ran with 2 threads (snakemake caps rule threads to the driver's
# CPU count when --cores is unset) so fastqc/trim_galore couldn't finish the
# big (1.3-2.1 GB) inputs within the 2-day genomics limit. Cancelled 2026-06-01.
#
# This resume run:
#   * NO --forceall  -> the 182 finished cells are skipped.
#   * --rerun-triggers mtime + --rerun-incomplete -> only the 6 missing cells
#     (85 jobs: demux -> mapping[Bhmem] -> bamprocess -> qc -> bisqc ->
#      bistools[BisSNP] -> methylation x8).
#   * --cores 64 raises the thread cap above the driver's 2 CPUs.
#   * --set-threads gives the 2 heaviest rules 16 threads (fastqc -t / trim -j),
#     which the slurm executor maps to --cpus-per-task=16. With 16 threads the
#     big cells finish in a few hours, well inside the 2-day genomics limit.
#   * --set-resources gives them LARGE memory (partition/runtime stay at the
#     profile default: genomics, 2 days).
#       demultiplex_fastqc_trim: 16 cpu, 64 GB, genomics (2-day default)
#       mapping                : 16 cpu, 48 GB, genomics (2-day default)
#
# Tooling (unchanged, from pipeline configs/config.yaml):
#   alignment   -> software/bisulfiteHic  (Bhmem, built 2026-05-11)
#   methylation -> Bis-SNP.latest.jar = BisSNP-1.2.bisqc_walkers.jar (2026-05-13)
#
# Submit:  sbatch codes/run_pipeline_gm_sc_new_resume.sh
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=scnh_resume
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnomehic/logs/pipeline_rerun/resume.%j.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnomehic/logs/pipeline_rerun/resume.%j.err

set -eo pipefail

source /projects/b1198/epifluidlab/yoshii/software/conda/etc/profile.d/conda.sh
conda activate scnomehic
module load java/jdk-17.0.2+8

PIPELINE=/gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline
cd "${PIPELINE}"

echo "=== Resume driver start: $(date) on $(hostname) ==="
which sbatch && sbatch --version

snakemake \
  -s "${PIPELINE}/Snakefile" \
  --profile "${PIPELINE}/profiles/slurm" \
  --configfile "${PIPELINE}/configs/config.yaml" \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --cores 64 \
  --set-threads \
      demultiplex_fastqc_trim=16 \
      mapping=16 \
  --set-resources \
      demultiplex_fastqc_trim:mem_mb=64000 \
      mapping:mem_mb=48000 \
  -j 1000 \
  -p

echo "=== Resume driver done: $(date) ==="
