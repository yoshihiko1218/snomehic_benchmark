#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=nomeqc
#SBATCH --output=logs/05.nome_qc/nomeqc.%a.out
#SBATCH --error=logs/05.nome_qc/nomeqc.%a.err

# Per-cell NOMe QC: detected HCG (CpG) and GCH (GpC) site counts + chrM/chr21
# trinucleotide methylation rates (ACT=conversion control, ACG=endogenous CpG,
# GCT=GpC accessibility), computed from the Bismark NOMe cov files and the
# deduplicated BAM (codes/nome_qc_sites_trinuc.py).
#
# LISTFILE selects which cells to process (one prefix per line); prefix must match
# the <prefix>_<mate>.rmdup.RG.bam and 05.methy/<prefix>_<mate>.NOMe.* naming.
#   GM   : sbatch --array=1-12 --export=ALL,LISTFILE=acc_list.txt            codes/04.qc_nome_sites.sh
#   K562 : sbatch --array=1-11 --export=ALL,LISTFILE=acc_list_k562_cellids.txt codes/04.qc_nome_sites.sh

source /home/jmj7858/.bashrc
conda activate scnomehic
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome || { echo "ERROR: cd failed"; exit 1; }
mkdir -p logs/05.nome_qc qc_stats

LISTFILE=${LISTFILE:-acc_list.txt}
REF=/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

prefix=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n' ${LISTFILE})
echo "[$(date)] cell=${prefix} (list=${LISTFILE})"
if [ -z "${prefix}" ]; then echo "ERROR: empty prefix for task ${SLURM_ARRAY_TASK_ID}"; exit 1; fi

python codes/nome_qc_sites_trinuc.py \
    --cell "${prefix}" \
    --methy_dir 05.methy \
    --align_dir 04.alignment \
    --ref "${REF}" \
    --chroms chrM,chr21 \
    --bam_suffix .rmdup.RG.bam \
    --out "qc_stats/${prefix}.nome_qc.tsv"
echo "[$(date)] DONE cell=${prefix}"
