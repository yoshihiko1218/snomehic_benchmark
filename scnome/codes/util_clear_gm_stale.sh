#!/bin/bash
# Remove STALE GM12878 outputs (SRR3729642-SRR3729653) that were produced with
# the OLD 5'-only trimming, so the re-trim / re-align / re-extract / re-QC
# pipeline regenerates them cleanly. K562 (SRR3729661-3729682) and control
# (SRR3729654-3729660) files are NOT touched.
#
# Usage:
#   bash codes/util_clear_gm_stale.sh         # DRY RUN: only lists what would be deleted
#   bash codes/util_clear_gm_stale.sh --yes   # actually delete
set -u

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

# GM12878 cells = first 12 lines of acc_list.txt (SRR3729642-SRR3729653)
GM=$(awk 'NR>=1 && NR<=12' acc_list.txt)
DIRS="03.trimmed_fastq 04.alignment 05.methy qc_stats"

DELETE=0
if [ "${1:-}" = "--yes" ]; then
  DELETE=1
fi

echo "GM cells:"
echo "${GM}" | tr '\n' ' '
echo ""
echo "Directories: ${DIRS}"
echo ""

n=0
for srr in ${GM}; do
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
