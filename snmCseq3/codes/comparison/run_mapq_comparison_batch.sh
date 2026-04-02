#!/bin/bash
# Run MAPQ comparison for all samples present in both bhmem and yap.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Samples with both bhmem BAM and yap 3C BAM (Group22, 42, 20, 34)
SAMPLES="SRR21549289 SRR21549292 SRR21549298 SRR21549299 SRR21549291"
for sample in ${SAMPLES}; do
  "${SCRIPT_DIR}/run_mapq_comparison.sh" "${sample}" || true
done
echo "Batch done. Check mapq_comparison/ under snmCseq3/"
