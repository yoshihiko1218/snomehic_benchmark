#!/bin/bash
# Aggregate per-cell HCG/GCH/HCH counts + methylation into one table.
# Run after the 05 array finishes. Per-cell file columns (tab):
#   cell HCGn HCGmc HCGcov GCHn GCHmc GCHcov HCHn HCHmc HCHcov
set -euo pipefail
BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
COUNT_DIR="${BASE}/mapping/nome_counts"
OUT="${BASE}/mapping/stats/hcg_gch_site_counts.tsv"

{
  printf "cell_id\tHCG_site_count\tGCH_site_count\tHCH_site_count\tHCG_mc_rate\tGCH_mc_rate\tHCH_mc_rate\n"
  cat "${COUNT_DIR}"/*.txt | sort | awk -F'\t' '{
    hcg_r = ($4>0)? 100*$3/$4 : 0;
    gch_r = ($7>0)? 100*$6/$7 : 0;
    hch_r = ($10>0)? 100*$9/$10 : 0;
    printf "%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n", $1, $2, $5, $8, hcg_r, gch_r, hch_r
  }'
} > "${OUT}"

n=$(($(wc -l < "${OUT}") - 1))
echo "Wrote ${OUT} with ${n} cells"
echo "=== medians ==="
tail -n +2 "${OUT}" | awk -F'\t' '{c[NR]=$2;g[NR]=$3; hr[NR]=$5; gr[NR]=$6; br[NR]=$7}
  END{
    asort(c); asort(g); asort(hr); asort(gr); asort(br); m=int(NR/2);
    printf "HCG_site_count median ~ %d ; GCH_site_count median ~ %d\n", c[m], g[m];
    printf "HCG_mc_rate median ~ %.1f%% ; GCH_mc_rate median ~ %.2f%% ; HCH(bg) median ~ %.2f%%\n", hr[m], gr[m], br[m];
  }'
