#!/bin/bash
# Compute Bismark-native trinuc QC on a few Smallwood ESC cells and compare to
# the existing BisSNP-derived values. Uses the SAME shared script as scNOMe
# (scnome/codes/nome_qc_sites_trinuc.py) with --mates "" (single no-mate cell).
source /home/jmj7858/.bashrc
conda activate scnomehic
set -o pipefail

cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood
QC=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome/codes/nome_qc_sites_trinuc.py
REF=/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa
OUT=qc_stats/trinuc_bismark
mkdir -p ${OUT}

for cell in SRR1248457 SRR1248458 SRR1248481 SRR1248490 SRR1248496; do
    echo "[`date`] === ${cell} ==="
    python3 ${QC} \
        --cell ${cell} \
        --methy_dir 06.methy/hcg \
        --align_dir 05.align_mm10 \
        --ref ${REF} \
        --chroms chrM,chr19 \
        --bam_suffix .rmdup.RG.bam \
        --mates "" \
        --out ${OUT}/${cell}.trinuc_bismark.csv
done
echo "[`date`] ALL DONE"
