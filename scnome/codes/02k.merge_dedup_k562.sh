#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-11
#SBATCH --job-name=k562merge
#SBATCH --output=logs/04.k562_merge/merge.%a.txt
#SBATCH --error=logs/04.k562_merge/merge.%a.err

# K562 cells were each sequenced as TWO runs (consecutive SRR pair = one library).
# To call methylation per cell we MERGE the two runs' bismark alignments and then
# remove duplicates ONCE on the merged data -- this removes PCR duplicates that
# occur WITHIN a run AND ACROSS the two runs (a fragment sequenced in both runs).
# Deduping each run separately would leave cross-run duplicates -> double counting.
#
# Reuses the existing per-run bismark BAMs (04.alignment/<run>.<run>_<m>_trimmed_bismark_bt2.bam),
# which were aligned with the correct K562 settings (5' clip, --non_directional,
# --score_min L,0,-0.2). markdup is run WITHOUT -r (marks, not removes), identical
# to how every other cell type was deduped -- only the trimming differs by type.

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome || { echo "ERROR: cd failed"; exit 1; }
conda activate scnomehic
set -o pipefail
mkdir -p logs/04.k562_merge

MAP=acc_list_k562_cells.tsv
ALIGN=04.alignment

read cell r1 r2 < <(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n{print $1, $2, $3}' ${MAP})
echo "[$(date)] cell=${cell} run1=${r1} run2=${r2}"
if [ -z "${cell}" ] || [ -z "${r1}" ] || [ -z "${r2}" ]; then echo "ERROR: bad map line ${SLURM_ARRAY_TASK_ID}"; exit 1; fi

RGSM="${cell}"; RGPL="ILLUMINA"; RGLB="lib1"; RGPU="unit1"

for m in 1 2; do
  b1=${ALIGN}/${r1}.${r1}_${m}_trimmed_bismark_bt2.bam
  b2=${ALIGN}/${r2}.${r2}_${m}_trimmed_bismark_bt2.bam
  for b in "${b1}" "${b2}"; do
    if [ ! -s "${b}" ]; then echo "ERROR: missing input bismark BAM ${b}"; exit 1; fi
  done

  echo "[run ] mate ${m}: cat ${r1}+${r2} | sort | markdup -> ${cell}_${m}.rmdup.bam"
  samtools cat "${b1}" "${b2}" \
    | samtools sort -@ 16 -T tmp_${cell}_${m} -O bam - \
    | samtools markdup -@ 16 - ${ALIGN}/${cell}_${m}.rmdup.bam
  if [ ! -s "${ALIGN}/${cell}_${m}.rmdup.bam" ]; then echo "ERROR: no merged rmdup.bam for mate ${m}"; exit 1; fi

  samtools addreplacerg \
    -r "@RG\tID:${cell}_${m}\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
    -o ${ALIGN}/${cell}_${m}.rmdup.RG.bam \
    ${ALIGN}/${cell}_${m}.rmdup.bam
  samtools index ${ALIGN}/${cell}_${m}.rmdup.RG.bam

  python /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/bam_summary_universal.py \
    --in_bam ${ALIGN}/${cell}_${m}.rmdup.RG.bam --out_summary ${ALIGN}/${cell}_${m}.summary.txt
done
echo "[$(date)] DONE cell=${cell}"
