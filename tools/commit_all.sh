#!/bin/bash
set -euo pipefail
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark

echo "=== Staging .gitignore ==="
git add .gitignore

echo "=== Staging top-level files ==="
git add PROJECT_CONTEXT.md 2>/dev/null || true
git add bam_summary_universal.py 2>/dev/null || true

echo "=== Staging snmCseq3 code files ==="
# codes directory - scripts and configs
find snmCseq3/codes -type f \( -name "*.py" -o -name "*.sh" -o -name "*.md" -o -name "*.ini" -o -name "*.yaml" -o -name "*.yml" \) ! -path "*__pycache__*" -exec git add {} \;

# cigar_nm_walk subdir
find snmCseq3/codes/comparison/cigar_nm_walk -type f \( -name "*.py" -o -name "*.sh" -o -name "*.md" \) 2>/dev/null -exec git add {} \; || true

# acc_list
git add snmCseq3/acc_list.txt

echo "=== Staging Snakefiles ==="
find snmCseq3/alignment -name "Snakefile" -exec git add {} \;
find snmCseq3/alignment_bowtie1 -name "Snakefile" -exec git add {} \;
find snmCseq3/alignment_mapq0 -name "Snakefile" -exec git add {} \;

# mapq0 config and stats
git add snmCseq3/alignment_mapq0/mapping_config.ini 2>/dev/null || true
find snmCseq3/alignment_mapq0/snakemake -type f -exec git add {} \; 2>/dev/null || true

# MappingSummary files (small compressed)
find snmCseq3/alignment -name "MappingSummary.csv.gz" -exec git add {} \; 2>/dev/null || true
find snmCseq3/alignment_mapq0 -name "MappingSummary.csv.gz" -exec git add {} \; 2>/dev/null || true

echo "=== Staging other benchmark dirs (code only) ==="
for dir in droplethic methylhic methylhic_new nagano scnome scnomehic smallwood snmCseq2 summary figures; do
    if [ -d "$dir" ]; then
        find "$dir" -type f \( -name "*.py" -o -name "*.sh" -o -name "*.R" -o -name "*.r" -o -name "*.md" -o -name "*.ini" -o -name "*.yaml" -o -name "*.yml" -o -name "*.toml" -o -name "*.json" -o -name "*.cfg" -o -name "Snakefile" -o -name "Makefile" \) ! -path "*__pycache__*" ! -path "*.snakemake*" ! -path "*/.cursor/*" -exec git add {} \; 2>/dev/null || true
    fi
done

# figures - add PDFs/PNGs/SVGs under 5MB
find figures -type f \( -name "*.pdf" -o -name "*.png" -o -name "*.svg" \) -size -5M -exec git add {} \; 2>/dev/null || true

echo "=== Summary ==="
git diff --cached --stat | tail -5
echo "Total files staged: $(git diff --cached --name-only | wc -l)"

echo "=== Committing ==="
git commit -m "Add all benchmark code, scripts, Snakefiles, and configs" \
  -m "Includes pipeline scripts, comparison analysis code, Snakefiles for" \
  -m "snmCseq3 (bowtie2, bowtie1, mapq0), and code for other benchmark datasets." \
  -m "" \
  -m "Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>"

echo "=== Pushing ==="
GIT_ASKPASS="" git push origin main

echo "=== Done ==="
