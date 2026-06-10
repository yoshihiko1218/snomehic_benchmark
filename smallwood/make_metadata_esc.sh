#!/usr/bin/env bash
# Build metadata_esc.tsv and subset the QC summary to the 32 ESC cells.
# Cell-type ranges provided by user (Smallwood 2014, GSE56879):
#   2i ESCs : SRR1248457 .. SRR1248468  -> 2i_1  .. 2i_12
#   Ser ESCs: SRR1248477 .. SRR1248496  -> Ser_1 .. Ser_20
set -euo pipefail
cd "$(dirname "$0")"

# --- metadata_esc.tsv ---
{
  printf 'SRR\tcondition\tlabel\n'
  i=1
  for srr in SRR12484{57..68}; do printf '%s\t2i\t2i_%d\n' "$srr" "$i"; i=$((i+1)); done
  i=1
  for srr in SRR12484{77..96}; do printf '%s\tSer\tSer_%d\n' "$srr" "$i"; i=$((i+1)); done
} > metadata_esc.tsv
echo "metadata_esc.tsv: $(( $(wc -l < metadata_esc.tsv) - 1 )) cells"

# --- subset QC summary (keep header + rows whose CellID is in acc_list_esc.txt) ---
awk -F, 'NR==FNR{keep[$1]=1; next} FNR==1 || ($1 in keep)' \
    acc_list_esc.txt smallwood_qc_summary.csv > smallwood_qc_summary_esc.csv
echo "smallwood_qc_summary_esc.csv: $(( $(wc -l < smallwood_qc_summary_esc.csv) - 1 )) data rows"
