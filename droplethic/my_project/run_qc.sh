#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1
#SBATCH --job-name=qc.%a
#SBATCH --output=logs/qc/qc.%a.out
#SBATCH --error=logs/qc/qc.%a.err

# ── EDIT THESE ──
# Set --array=1-N where N = number of samples in config.yaml
# e.g. for 5 samples: --array=1-5

source ~/.bashrc
conda activate rupture

PROJECT_DIR=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/droplethic/my_project
GENOME=hg38
MAPQ=30

cd ${PROJECT_DIR}
mkdir -p logs/qc qc_stats

# Get sample ID from config.yaml using the array task ID
# Extracts the Nth sample_id entry (skipping the "sample_id:" line)
SAMPLE=$(grep -A 1000 'sample_id:' config.yaml \
    | grep '^\s*-' \
    | sed 's/.*-\s*//' \
    | awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}')

echo "=== Processing sample: ${SAMPLE} ==="
echo "BAM: 03.mapping/${SAMPLE}_${GENOME}.bam"
echo "Pairs: 03.mapping/${SAMPLE}_${GENOME}.sc.pairs.gz"

# Build arguments based on what files exist
ARGS="--output qc_stats/${SAMPLE}_${GENOME} --mapq ${MAPQ}"

if [ -f "03.mapping/${SAMPLE}_${GENOME}.bam" ]; then
    ARGS="${ARGS} --bam 03.mapping/${SAMPLE}_${GENOME}.bam"
else
    echo "WARNING: BAM not found, skipping BAM-based stats"
fi

if [ -f "03.mapping/${SAMPLE}_${GENOME}.sc.pairs.gz" ]; then
    ARGS="${ARGS} --pairs 03.mapping/${SAMPLE}_${GENOME}.sc.pairs.gz"
else
    echo "WARNING: Pairs file not found, skipping pairs-based stats"
fi

python rupture_qc.py ${ARGS}

echo "=== Done: ${SAMPLE} ==="
