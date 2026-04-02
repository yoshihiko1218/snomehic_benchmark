#!/bin/bash
# Run MAPQ comparison for samples present in both bhmem and yap pipelines.
# Usage: ./run_mapq_comparison.sh [sample_id]
#   If sample_id omitted, runs on SRR21549292 (in Group22).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BHMEM_BAM="${BASE}/04.bhmem_bam"
YAP_ALIGN="${BASE}/alignment"
OUT="${BASE}/mapq_comparison"
ACC_LIST="${BASE}/acc_list.txt"

mkdir -p "${OUT}"

# Map sample -> yap Group (from grep of CELL_IDS)
find_group() {
  local sid="$1"
  for g in "${YAP_ALIGN}"/Group*/Snakefile; do
    if grep -q "'${sid}'" "${g}" 2>/dev/null; then
      echo "$(dirname "${g}")"
      return
    fi
  done
  echo ""
}

run_sample() {
  local prefix="$1"
  local bhmem_bam="${BHMEM_BAM}/${prefix}.bhmem.bam"
  local group_dir
  group_dir="$(find_group "${prefix}")"
  if [[ -z "${group_dir}" ]]; then
    echo "SKIP: ${prefix} not in any yap Group" >&2
    return 1
  fi
  local yap_bam="${group_dir}/bam/${prefix}.3C.sorted.bam"
  if [[ ! -f "${bhmem_bam}" ]]; then
    echo "SKIP: bhmem BAM missing: ${bhmem_bam}" >&2
    return 1
  fi
  if [[ ! -f "${yap_bam}" ]]; then
    echo "SKIP: yap 3C BAM missing: ${yap_bam}" >&2
    return 1
  fi

  echo "Processing: ${prefix}"
  local sample_out="${OUT}/${prefix}"
  mkdir -p "${sample_out}"

  python3 "${SCRIPT_DIR}/extract_mapq.py" bhmem "${bhmem_bam}" -o "${sample_out}/bhmem_mapq.tsv"
  python3 "${SCRIPT_DIR}/extract_mapq.py" yap   "${yap_bam}"   -o "${sample_out}/yap_mapq.tsv" --primary-only

  python3 "${SCRIPT_DIR}/compare_mapq.py" \
    "${sample_out}/bhmem_mapq.tsv" \
    "${sample_out}/yap_mapq.tsv" \
    -o "${sample_out}/mapq_comparison" \
    --plot

  # Subset: yap>30, bhmem<30
  awk -F'\t' 'NR>1 && $3<30 && $4>30' "${sample_out}/mapq_comparison.joined.tsv" | \
    (echo -e "base_id\tis_r1\tmapq_bhmem\tmapq_yap\tdiff"; cat) > "${sample_out}/yap_high_bhmem_low.tsv"

  # Summarize discrepant reads (chr, pos, distance, NM, etc.)
  python3 "${SCRIPT_DIR}/summarize_discrepant_reads.py" \
    "${sample_out}/yap_high_bhmem_low.tsv" \
    "${bhmem_bam}" \
    "${yap_bam}" \
    -o "${sample_out}/discrepant_summary" \
    --mapq-tsv "${sample_out}/yap_high_bhmem_low.tsv"

  echo "Done: ${prefix} -> ${sample_out}/"
}

if [[ $# -ge 1 ]]; then
  run_sample "$1"
else
  # Default: SRR21549292 (in Group22, has bhmem BAM)
  run_sample "SRR21549292"
fi
