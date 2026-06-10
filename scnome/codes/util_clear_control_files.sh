#!/bin/bash
# Delete ALL control / spike-in sample files (SRR3729654-SRR3729660). These
# 7 samples are excluded from analysis and have been removed from acc_list.txt
# (original preserved as acc_list_with_controls.txt). The control SRR IDs are
# hardcoded here because they are no longer present in acc_list.txt.
#
# Usage:
#   bash codes/util_clear_control_files.sh         # DRY RUN: only lists what would be deleted
#   bash codes/util_clear_control_files.sh --yes   # actually delete
set -u

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

CONTROLS="SRR3729654 SRR3729655 SRR3729656 SRR3729657 SRR3729658 SRR3729659 SRR3729660"
DIRS="01.fastq 02.fastqc_out 03.trimmed_fastq 04.alignment 05.methy qc_stats"

DELETE=0
if [ "${1:-}" = "--yes" ]; then
  DELETE=1
fi

echo "Control cells: ${CONTROLS}"
echo "Directories:   ${DIRS}"
echo ""

n=0
for srr in ${CONTROLS}; do
  for d in ${DIRS}; do
    for f in ${d}/${srr}*; do
      [ -e "${f}" ] || continue
      n=$((n + 1))
      if [ "${DELETE}" -eq 1 ]; then
        rm -f "${f}"
        echo "DELETED ${f}"
      else
        echo "would delete  ${f}"
      fi
    done
  done
done

echo ""
echo "Total matched files: ${n}"
if [ "${DELETE}" -eq 1 ]; then
  echo "Deletion complete."
else
  echo "DRY RUN only. Re-run with --yes to delete."
fi
