#!/bin/bash
# Isolate WHY Bismark-native != BisSNP: re-run SRR1248481 chr19 after applying
# BisSNP-like filters (MAPQ>=30 + exclude duplicates). If endo jumps toward the
# BisSNP value (77.606), filtering is the driver (not a context/strand bug).
source /home/jmj7858/.bashrc
conda activate scnomehic
set -o pipefail
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/smallwood
SAM=/gpfs/projects/b1198/epifluidlab/yoshii/software/samtools-1.16/bin/samtools
QC=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome/codes/nome_qc_sites_trinuc.py
REF=/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa
cell=SRR1248481
mkdir -p qc_stats/trinuc_filter_test

# MAPQ>=30 + exclude dups (0x400), like BisSNP --mmq 30 + dup filter
${SAM} view -b -q 30 -F 1024 -o qc_stats/trinuc_filter_test/${cell}.q30.nodup.bam 05.align_mm10/${cell}.rmdup.RG.bam
${SAM} index qc_stats/trinuc_filter_test/${cell}.q30.nodup.bam

python3 ${QC} \
    --cell ${cell} \
    --methy_dir 06.methy/hcg \
    --align_dir qc_stats/trinuc_filter_test \
    --ref ${REF} \
    --chroms chr19 \
    --bam_suffix .q30.nodup.bam \
    --mates "" \
    --out qc_stats/trinuc_filter_test/${cell}.filtered.csv
echo "=== BisSNP ref: chr19_endo=77.606 chr19_noncpg=0.648 ==="
echo "=== raw Bismark (earlier): chr19_endo=57.048 chr19_noncpg=2.752 ==="
echo "ALL DONE"
