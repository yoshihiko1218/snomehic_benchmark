#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100
#SBATCH --job-name=nome_hcg_gch
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/05.nome/nome.%a.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/05.nome/nome.%a.err

# Recompute ALLC in NOMe context (num_upstr_bases=1) from the retained deduped DNA
# BAMs, then count covered HCG and GCH loci (cov>=1) per cell.
#   HCG = H-CG (CpG not preceded by G): contexts ACG/CCG/TCG
#   GCH = G-CH (GpC not followed by G): contexts GCA/GCC/GCT
#   GCG excluded from both (standard NOMe convention).

# source bashrc / activate conda BEFORE set -u (bashrc has an unbound var)
source /home/jmj7858/.bashrc
conda activate mapping
export PATH="$CONDA_PREFIX/bin:$PATH"
set -euo pipefail

BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
cd "${BASE}"   # ensure a stable cwd; manifest also uses absolute paths
REF_FASTA="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
MANIFEST="${BASE}/codes/nome_manifest.tsv"
COUNT_DIR="${BASE}/mapping/nome_counts"
mkdir -p "${COUNT_DIR}" "${BASE}/codes/logs/05.nome"

# Nth manifest line: cell_id <tab> bam_path <tab> allc_out
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")
cell=$(echo "${line}" | cut -f1)
bam=$(echo "${line}" | cut -f2)
allc_out=$(echo "${line}" | cut -f3)

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: cell=${cell}"
echo "bam=${bam}"; echo "allc_out=${allc_out}"

# 1) NOMe-context ALLC (skip if already present and non-empty)
if [[ ! -s "${allc_out}" ]]; then
    allcools bam-to-allc \
        --bam_path "${bam}" \
        --reference_fasta "${REF_FASTA}" \
        --output_path "${allc_out}" \
        --num_upstr_bases 1 \
        --num_downstr_bases 1 \
        --cpu 1 \
        --compress_level 5
fi

# 2) Per-context covered loci + methylation (mc/cov). Column 4 = trinucleotide context,
#    col5 = methylated count (mc), col6 = total coverage (cov).
#    HCG = ACG/CCG/TCG ; GCH = GCA/GCC/GCT ; HCH = H C H (background, excludes any CpG/GpC).
zcat "${allc_out}" | awk -F'\t' -v cell="${cell}" '
    function isH(b){ return (b=="A"||b=="C"||b=="T") }
    {
        c=$4; mc=$5; cov=$6;
        u=substr(c,1,1); d=substr(c,3,1);
        if ((u=="A"||u=="C"||u=="T") && d=="G") { hcg_n++; hcg_mc+=mc; hcg_cov+=cov }
        else if (u=="G" && (d=="A"||d=="C"||d=="T")) { gch_n++; gch_mc+=mc; gch_cov+=cov }
        else if ((u=="A"||u=="C"||u=="T") && (d=="A"||d=="C"||d=="T")) { hch_n++; hch_mc+=mc; hch_cov+=cov }
    }
    END {
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            cell, hcg_n, hcg_mc, hcg_cov, gch_n, gch_mc, gch_cov, hch_n, hch_mc, hch_cov
    }' > "${COUNT_DIR}/${cell}.txt"

echo "[$(date)] DONE ${cell}:"; cat "${COUNT_DIR}/${cell}.txt"
