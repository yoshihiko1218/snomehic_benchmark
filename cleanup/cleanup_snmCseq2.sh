#!/usr/bin/env bash
# Cleanup for snmCseq2 — REVIEW BEFORE RUNNING. Run by hand from the benchmark root.
# KEEPS both run arms: Bismark-SE (05.align + 06.methy cov/reports -> snmcseq2_qc_summary.csv)
# AND the user's YAP run (yap_mapping_hg38, yap_mapping_mm10). Deletes only regenerable
# bismark_methylation_extractor context-txt intermediates + small abandoned dirs.
# Verified 2026-06-29: context txt unreferenced; cov.gz/trinuc/qc_summary already exist.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1   # -> benchmark/
B=snmCseq2

# --- BIG: per-read methylation context dumps (~612G), distilled output already exists ---
rm -f "$B"/06.methy/CHH_context_*.txt     # 447G
rm -f "$B"/06.methy/CHG_context_*.txt     # 136G
rm -f "$B"/06.methy/CpG_context_*.txt     # 30G

# --- small dead-ends ---
rm -rf "$B/.old_single_genome_attempt"    # 1.7M abandoned mm10-only YAP run

# --- OPTIONAL (uncomment if you want): regenerable post-trim FastQC reports ---
# rm -rf "$B/04.fastqc_out_2"             # 313M

echo "snmCseq2 cleanup done."
