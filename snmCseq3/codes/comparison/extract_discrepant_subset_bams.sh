#!/bin/bash
# Extract the same read base_ids (first N from yap_high_bhmem_low.tsv) into
# sorted, indexed BAMs for bhmem, yap Bowtie2, and yap Bowtie1.
#
# Usage:
#   ./extract_discrepant_subset_bams.sh [N]
#   N defaults to 10. Override BASE if needed (snmCseq3 root).
#
# Output (under mapq_comparison/SRR21549292/subset_bams/):
#   discrepant_subset.bhmem.sorted.bam(.bai)
#   discrepant_subset.yap_bowtie2.sorted.bam(.bai)
#   discrepant_subset.yap_bowtie1.sorted.bam(.bai)
#   discrepant_subset_read_ids.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"
N="${1:-10}"

SAMPLE="SRR21549292"
TSV="${BASE}/mapq_comparison/${SAMPLE}/yap_high_bhmem_low.tsv"
OUTDIR="${BASE}/mapq_comparison/${SAMPLE}/subset_bams"
READLIST="${OUTDIR}/discrepant_subset_read_ids.txt"

BH="${BASE}/04.bhmem_bam/${SAMPLE}.bhmem.bam"
BT2="${BASE}/alignment/Group22/bam/${SAMPLE}.3C.sorted.bam"
BT1="${BASE}/alignment_bowtie1/Group25/bam/${SAMPLE}.3C.sorted.bam"

for f in "$TSV" "$BH" "$BT2"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing: $f" >&2
    exit 1
  fi
done
if [[ ! -f "$BT1" ]]; then
  echo "WARN: Bowtie1 BAM not found (use GPFS copy if needed): $BT1" >&2
  echo "      Skipping yap_bowtie1 subset." >&2
  SKIP_BT1=1
else
  SKIP_BT1=0
fi

mkdir -p "$(dirname "$READLIST")"
mkdir -p "$OUTDIR"

awk -v n="$N" \
  'NR>1 && !seen[$1]++ { print $1; if (++c>=n) exit }' "$TSV" > "$READLIST"

echo "Read IDs ($(wc -l < "$READLIST") lines) -> $READLIST"

# --- bhmem: exact QNAME match ---
samtools view -h -N "$READLIST" "$BH" \
  | samtools sort -o "${OUTDIR}/discrepant_subset.bhmem.sorted.bam" -
samtools index "${OUTDIR}/discrepant_subset.bhmem.sorted.bam"
echo "  OK: discrepant_subset.bhmem.sorted.bam"

# --- yap (Bowtie2 & Bowtie1): QNAME starts with base_id + "_" ---
awk_filter() {
  # stdin: SAM lines; READLIST file in awk
  awk -v f="$READLIST" '
    BEGIN { while ((getline l < f) > 0) n[l]=1 }
    /^@/ { print; next }
    {
      for (k in n) {
        if (index($1, k "_") == 1) { print; next }
      }
    }
  '
}

samtools view -h "$BT2" | awk_filter \
  | samtools sort -o "${OUTDIR}/discrepant_subset.yap_bowtie2.sorted.bam" -
samtools index "${OUTDIR}/discrepant_subset.yap_bowtie2.sorted.bam"
echo "  OK: discrepant_subset.yap_bowtie2.sorted.bam"

if [[ "$SKIP_BT1" -eq 0 ]]; then
  samtools view -h "$BT1" | awk_filter \
    | samtools sort -o "${OUTDIR}/discrepant_subset.yap_bowtie1.sorted.bam" -
  samtools index "${OUTDIR}/discrepant_subset.yap_bowtie1.sorted.bam"
  echo "  OK: discrepant_subset.yap_bowtie1.sorted.bam"
fi

echo "Done. Output in: $OUTDIR"
