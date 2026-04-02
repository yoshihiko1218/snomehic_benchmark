#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1379
#SBATCH --job-name=qc_trim
#SBATCH --output=logs/01.qc_trim/qc_trim.%a.txt
#SBATCH --error=logs/01.qc_trim/qc_trim.%a.err

# set -euo pipefail

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3/

conda activate scnomehic

ACC_LIST="acc_list.txt"
INFOLDER="fastq"
OUTFOLDER_QC="02.fastqc_out"
OUTFOLDER_FQ="03.trimmed_fastq"

mkdir -p "${OUTFOLDER_QC}" "${OUTFOLDER_FQ}" logs/01.qc_trim

OVERLAP=6
QTHRESH=20
LTHRESH=30

R1_LEFT_CUT=10
R1_RIGHT_CUT=10
R2_LEFT_CUT=10
R2_RIGHT_CUT=10
R1_ADAPTER="AGATCGGAAGAGCACACGTCTGAAC"
R2_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGA"

THREADS="${SLURM_CPUS_PER_TASK}"

# --------------------
# Pick sample for this array task
# --------------------
prefix=$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' "${ACC_LIST}")
if [[ -z "${prefix}" ]]; then
  echo "ERROR: Could not find prefix for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${ACC_LIST}" >&2
  exit 1
fi

R1_IN="${INFOLDER}/${prefix}-R1.fq.gz"
R2_IN="${INFOLDER}/${prefix}-R2.fq.gz"

if [[ ! -s "${R1_IN}" || ! -s "${R2_IN}" ]]; then
  echo "ERROR: Missing input fastq: ${R1_IN} or ${R2_IN}" >&2
  exit 1
fi

echo "Sample: ${prefix}"
echo "R1: ${R1_IN}"
echo "R2: ${R2_IN}"

# --------------------
# FastQC on raw reads
# --------------------
fastqc -t "${THREADS}" --outdir "${OUTFOLDER_QC}" "${R1_IN}" "${R2_IN}"
echo "Raw FastQC done."

# --------------------
# Trimming (match Snakefile cutadapt logic)
#   Pass 1: adapter trim (-a)
#   Pass 2: overlap + quality + fixed clip + min length
# --------------------

# Output files (Snakefile-like naming)
R1_OUT="${OUTFOLDER_FQ}/${prefix}-R1.trimmed.fq.gz"
R2_OUT="${OUTFOLDER_FQ}/${prefix}-R2.trimmed.fq.gz"

R1_STATS="${OUTFOLDER_FQ}/${prefix}-R1.trimmed.stats.tsv"
R2_STATS="${OUTFOLDER_FQ}/${prefix}-R2.trimmed.stats.tsv"

# safety: avoid silent overwrite
for f in "${R1_OUT}" "${R2_OUT}" "${R1_STATS}" "${R2_STATS}"; do
  if [[ -e "${f}" ]]; then
    echo "ERROR: Output exists (won't overwrite): ${f}" >&2
    exit 1
  fi
done

# R1
cutadapt -a "${R1_ADAPTER}" "${R1_IN}" 2> "${R1_STATS}" \
| cutadapt  -O "${OVERLAP}" -q "${QTHRESH}" \
    -u "${R1_LEFT_CUT}" -u "-${R1_RIGHT_CUT}" -m "${LTHRESH}" \
    -o "${R1_OUT}" - >> "${R1_STATS}"

# R2
cutadapt --report=minimal -a "${R2_ADAPTER}" "${R2_IN}" 2> "${R2_STATS}" \
| cutadapt --report=minimal -O "${OVERLAP}" -q "${QTHRESH}" \
    -u "${R2_LEFT_CUT}" -u "-${R2_RIGHT_CUT}" -m "${LTHRESH}" \
    -o "${R2_OUT}" - >> "${R2_STATS}"

echo "Trimming done."
echo "R1 trimmed: ${R1_OUT}"
echo "R2 trimmed: ${R2_OUT}"

# --------------------
# FastQC on trimmed reads (optional but mirrors your example)
# --------------------
fastqc -t "${THREADS}" --outdir "${OUTFOLDER_QC}" "${R1_OUT}" "${R2_OUT}"
echo "Trimmed FastQC done."
