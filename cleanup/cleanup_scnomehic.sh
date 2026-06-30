#!/usr/bin/env bash
# Cleanup for scnomehic — REVIEW BEFORE RUNNING. Run by hand from the benchmark root.
# Deletes ONLY the bowtie1 parameter-experiment alignment dir, which is NOT read by
# summary/qc.ipynb or any summary/ script (verified 2026-06-29 — 0 active references).
# KEEPS: alignment/ (YAP/bowtie2, canonical — stats/MappingSummary.csv.gz consumed),
#        fastq/ (raw input), codes/, logs/, acc_list.txt, JOBS.md, SESSION_NOTE.
# The valid bhmem (method) output is EXTERNAL on b1198 and is untouched by this script.
#
# Frees ~341 GB.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1   # -> benchmark/

rm -rf scnomehic/alignment_bowtie1   # 341G  bowtie1 vs bowtie2 parameter experiment

echo "scnomehic cleanup done."
