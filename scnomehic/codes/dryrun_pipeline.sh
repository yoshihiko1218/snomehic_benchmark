#!/bin/bash
# Dry-run of the updated sc_NOMeHiC_pipeline against gm_sc_new (hg38).
# Shows what the default rule all (trim -> align -> bamproc -> bisqc -> methylation,
# NO hicluster) would run, given existing outputs. Resume-aware via mtime triggers.
set -eo pipefail

source /projects/b1198/epifluidlab/yoshii/software/conda/etc/profile.d/conda.sh
conda activate scnomehic
module load java/jdk-17.0.2+8 2>/dev/null || true

PIPELINE=/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline
cd "${PIPELINE}"

snakemake \
  -s "${PIPELINE}/Snakefile" \
  --configfile "${PIPELINE}/configs/config.yaml" \
  --rerun-triggers mtime \
  -np 2>&1 | tail -80
