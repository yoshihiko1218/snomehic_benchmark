#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --job-name=frag
#SBATCH --output=summary/frag_jobs/logs/%x.%A_%a.out
#SBATCH --error=summary/frag_jobs/logs/%x.%A_%a.err
# Uniform fragment-level QC: markdup-if-needed -> frag_counts, per cell.
# Usage: sbatch --array=1-N frag_array.sh <DATASET>
source /home/jmj7858/.bashrc 2>/dev/null
conda activate scnomehic 2>/dev/null
set -uo pipefail
B=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
cd "$B" || exit 1

DS="$1"
manifest="$B/summary/frag_jobs/${DS}.manifest.tsv"
outdir="$B/summary/frag_counts/${DS}"
mkdir -p "$outdir" "$B/summary/frag_jobs/logs"

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")
cell=$(echo "$line" | cut -f1)
bams=$(echo "$line" | cut -f2)
needs=$(echo "$line" | cut -f3)

tmpd="${SLURM_TMPDIR:-/tmp}/frag.${SLURM_ARRAY_JOB_ID:-x}.${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpd"
trap 'rm -rf "$tmpd"' EXIT

use_bams=""
IFS=',' read -ra arr <<< "$bams"
i=0
for bam in "${arr[@]}"; do
  if [ "$needs" = "1" ]; then
    samtools sort -@ 4 -T "$tmpd/s$i" -o "$tmpd/$i.sorted.bam" "$bam"
    samtools markdup -@ 4 "$tmpd/$i.sorted.bam" "$tmpd/$i.markdup.bam"
    use_bams="$use_bams $tmpd/$i.markdup.bam"
  else
    use_bams="$use_bams $bam"
  fi
  i=$((i+1))
done

python "$B/summary/frag_counts.py" --cell "$cell" --bam $use_bams > "$outdir/${cell}.tsv"
echo "done $DS $cell"
