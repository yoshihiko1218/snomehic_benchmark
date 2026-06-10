#!/bin/bash
# Fix yap mc-mode bug: pipelines/mc.py drops `bismark_reference` from the Snakefile
# header whenever the config still contains a `hisat3n_dna_reference` key (even a
# placeholder). `yap default-mapping-config` leaves that placeholder line in, so the
# generated mc.Snakefile (which calls `bismark {bismark_reference}`) hits a NameError.
# Fix = delete the hisat3n_dna_reference line so yap keeps bismark_reference.
set -euo pipefail

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq2"

CONFIGS=(
  "${BASE}/codes/mapping_config_yap_hg38.ini"
  "${BASE}/codes/mapping_config_yap_mm10.ini"
  "${BASE}/yap_mapping_hg38/mapping_config.ini"
  "${BASE}/yap_mapping_mm10/mapping_config.ini"
)

for cfg in "${CONFIGS[@]}"; do
  # Remove the hisat3n_dna_reference assignment line and its trailing comment line.
  grep -v -E '^hisat3n_dna_reference\s*=' "${cfg}" \
    | grep -v -E '^; reference prefix for the HISAT-3N DNA mapping' > "${cfg}.tmp"
  mv "${cfg}.tmp" "${cfg}"
  echo "patched: ${cfg}  (hisat3n_dna_reference present? $(grep -c hisat3n_dna_reference "${cfg}"))"
done
