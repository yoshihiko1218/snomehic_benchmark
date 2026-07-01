#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --job-name=s3ctc30
#SBATCH --output=summary/frag_jobs/logs/%x.%A_%a.out
#SBATCH --error=summary/frag_jobs/logs/%x.%A_%a.err
# snmCseq3 contacts at MapQ>=30 (in addition to YAP default MapQ>=10):
# filter the YAP 3C.sorted.bam at q30 and re-run yap-internal generate-contacts (min_gap=1000).
# Emits CisShort,CisLong(>1kb),Trans per cell.
source /home/jmj7858/.bashrc 2>/dev/null
conda activate mapping 2>/dev/null
set -uo pipefail
B=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
cd "$B" || exit 1
CS=/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.chrom.sizes
outdir="$B/summary/frag_counts/snmCseq3_contacts_q30"
mkdir -p "$outdir" "$B/summary/frag_jobs/logs"

manifest="$B/summary/frag_jobs/snmCseq3_contacts.manifest.tsv"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")
cell=$(echo "$line" | cut -f1)
bam=$(echo "$line" | cut -f2)

tmpd="${SLURM_TMPDIR:-/tmp}/s3c.${SLURM_ARRAY_JOB_ID:-x}.${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpd"; trap 'rm -rf "$tmpd"' EXIT

samtools view -b -q 30 "$bam" > "$tmpd/q30.bam"
yap-internal generate-contacts --bam_path "$tmpd/q30.bam" \
  --output_path "$tmpd/q30.contact.tsv.gz" --chrom_size_path "$CS" --min_gap 1000
# counts file -> CisShortContact,CisLongContact,TransContact
cp "$tmpd/q30.contact.tsv.counts.txt" "$outdir/${cell}.q30.counts.txt"
echo "done $cell : $(tr '\n' ' ' < $outdir/${cell}.q30.counts.txt)"
