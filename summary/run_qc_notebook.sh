#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=run_qc_nb
#SBATCH --output=summary/logs/run_qc_nb.out
#SBATCH --error=summary/logs/run_qc_nb.err

# Execute qc.ipynb end-to-end to regenerate all figures with the corrected data.
# Run a COPY at the benchmark root so the kernel CWD == benchmark/ and the notebook's
# relative paths (droplethic/, summary/, figures/) resolve correctly.
source /home/jmj7858/.bashrc
conda activate scnomehic
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
mkdir -p summary/logs

cp summary/qc.ipynb _tmp_qc_run.ipynb
echo "[$(date)] executing notebook ..."
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=1800 \
    _tmp_qc_run.ipynb
rc=$?
echo "[$(date)] nbconvert exit=$rc"
exit $rc
