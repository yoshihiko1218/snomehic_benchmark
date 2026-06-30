#!/usr/bin/env bash
# Cleanup for snmCseq3 — REVIEW BEFORE RUNNING. Run by hand from the benchmark root.
# Deletes ONLY bhmem-vs-YAP methodology-investigation dead-ends that are NOT
# consumed by summary/qc.ipynb or any summary/ script (verified 2026-06-29).
# KEEPS the canonical run: fastq, 03.trimmed_fastq, alignment (YAP),
# 04.bhmem_bam (bhmem), 02.fastqc_out, codes, logs, and small conclusion dirs.
#
# Frees ~275 GB. Each rm is listed separately so you can run a subset.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1   # -> benchmark/
B=snmCseq3

# --- large alignment experiments ---
rm -rf "$B/alignment_mapq0"          # 69G  MAPQ0 experiment
rm -rf "$B/alignment_bowtie1"        # 49G  bowtie1 vs bowtie2 experiment

# --- bhmem JAR / param / subset experiments ---
rm -rf "$B/04.bhmem_bam_subset"      # 31G
rm -rf "$B/04.bhmem_bam_xg"          # 17G  proven identical to base
rm -rf "$B/04.bhmem_bam_original"    # 17G  proven identical to base
rm -rf "$B/04.bhmem_bam_oldjar"      # 232M
rm -rf "$B/04.bhmem_bam_buf100k"     # 232M
rm -rf "$B/04.bhmem_bam_noxg"        # 220M

# --- test/subset inputs ---
rm -rf "$B/test_fastq"               # 24G
rm -rf "$B/test_fastq_renamed"       # 21G
rm -rf "$B/03.trimmed_fastq_subset"  # 24G

# --- superseded comparison outputs (kept: _v2, final_, comprehensive_) ---
rm -rf "$B/mapq_comparison_mapq0"    # 1.5G
rm -rf "$B/corrected_nm_comparison"  # 348M  round-1, superseded by _v2
rm -rf "$B/mapq_comparison"          # 269M

# --- loose scratch test BAMs at snmCseq3 root ---
rm -f "$B"/one.bhmem.bam "$B"/one.yap3c.bam \
      "$B"/one.yap_bowtie1.sorted.bam "$B"/one.yap_bowtie1.sorted.bam.bai \
      "$B"/pdf_subset.bhmem.bam "$B"/pdf_subset.bhmem.sorted.bam "$B"/pdf_subset.bhmem.sorted.bam.bai \
      "$B"/pdf_subset.yap3c.bam "$B"/pdf_subset.yap3c.sorted.bam "$B"/pdf_subset.yap3c.sorted.bam.bai \
      "$B"/pdf_subset.yap_bowtie1.sorted.bam "$B"/pdf_subset.yap_bowtie1.sorted.bam.bai

echo "snmCseq3 cleanup done."
