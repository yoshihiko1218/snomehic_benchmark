#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1379
#SBATCH --job-name=rename_gzip_fastq
#SBATCH --output=logs/00.rename_gzip/rename_gzip.%a.txt
#SBATCH --error=logs/00.rename_gzip/rename_gzip.%a.err

# set -euo pipefail

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1198/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3/fastq

# If acc_list.txt lives one dir above fastq/, adjust path accordingly
ACC_LIST="../acc_list.txt"

prefix="$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print $1; exit}' "${ACC_LIST}")"
if [[ -z "${prefix}" ]]; then
  echo "ERROR: No prefix found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi
echo "Task ${SLURM_ARRAY_TASK_ID}: prefix=${prefix}"

in1="${prefix}_1.fastq"
in2="${prefix}_2.fastq"
out1="${prefix}-R1.fq.gz"
out2="${prefix}-R2.fq.gz"

# Sanity checks
# if [[ ! -s "${in1}" ]]; then
#   echo "ERROR: Missing input: ${in1}" >&2
#   exit 1
# fi
if [[ ! -s "${in2}" ]]; then
  echo "ERROR: Missing input: ${in2}" >&2
  exit 1
fi
# if [[ -e "${out1}" || -e "${out2}" ]]; then
#   echo "ERROR: Output exists, not overwriting: ${out1} or ${out2}" >&2
#   exit 1
# fi

# Optional: lock per-prefix to prevent accidental double-processing
lockdir=".lock_${prefix}"
if ! mkdir "${lockdir}" 2>/dev/null; then
  echo "ERROR: Lock exists (${lockdir}). Another task may be processing ${prefix}." >&2
  exit 1
fi
trap 'rmdir "${lockdir}" >/dev/null 2>&1 || true' EXIT

echo "Compressing -> ${out1}, ${out2}"
# gzip -c "${in1}" > "${out1}"
gzip -c "${in2}" > "${out2}"

# Verify gz files are readable (cheap integrity check)
# gzip -t "${out1}"
gzip -t "${out2}"

# Remove originals only after successful gzip + test
rm -f "${in1}" "${in2}"

echo "DONE: ${prefix}"
