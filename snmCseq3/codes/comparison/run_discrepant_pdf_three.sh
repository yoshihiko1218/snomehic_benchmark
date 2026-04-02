#!/bin/bash
# Three-way PDF: bhmem + yap Bowtie2 + yap Bowtie1 (SRR21549292).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Bowtie1 run is Group25 on GPFS; Bowtie2 is Group22 (same sample)
BT1="${BASE}/alignment_bowtie1/Group25/bam/SRR21549292.3C.sorted.bam"
BT2="${BASE}/alignment/Group22/bam/SRR21549292.3C.sorted.bam"
BH="${BASE}/04.bhmem_bam/SRR21549292.bhmem.bam"
TSV="${BASE}/mapq_comparison/SRR21549292/yap_high_bhmem_low.tsv"
OUT="${BASE}/mapq_comparison/SRR21549292/discrepant_reads_three.pdf"

for f in "$BT1" "$BT2" "$BH" "$TSV"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing: $f" >&2
    exit 1
  fi
done

conda run -n scnomehic python3 "${SCRIPT_DIR}/make_discrepant_pdf_three.py" \
  "$TSV" "$BH" "$BT2" "$BT1" -o "$OUT" -n 10

echo "Wrote $OUT"
