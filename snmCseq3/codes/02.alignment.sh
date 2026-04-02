#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-100
#SBATCH --job-name=bhmem
#SBATCH --output=logs/02.bhmem/bhmem.%A_%a.out
#SBATCH --error=logs/02.bhmem/bhmem.%A_%a.err

#set -euo pipefail

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1198/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3
module load java/jdk-17.0.2+8
conda activate scnomehic
# -----------------------
# Paths
# -----------------------
ACC_LIST="acc_list.txt"

# Trimmed FASTQ directory (from trimming step)
TRIMDIR="03.trimmed_fastq"

# Output directory for BAM
OUTDIR="04.bhmem_bam"
mkdir -p "${OUTDIR}" logs/02.bhmem

# Reference FASTA (Bhmem expects parent dir contains Bisulfite_Genome/*conversion indexes)
REF_FA="/gpfs/projects/b1198/epifluidlab/shared/data/genomes/mm10/mm10.fa"

# Restriction enzyme bedgraph
ENZYME_BED="/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/dpnII.span_region.bedgraph"

# Bhmem classpath / JNI library path (these match your working methylhic call)
JNI_LIB="/home/jmj7858/epifluidlab/software/bisulfitehic/jbwa/jbwa-1.0.0/src/main/native"
CP="/home/jmj7858/epifluidlab/software/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/home/jmj7858/epifluidlab/software/bisulfitehic/jbwa/jbwa-1.0.0/jbwa.jar"
MAIN_CLASS="main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem"

# -----------------------
# Threads / params
# -----------------------
THREADS="${SLURM_CPUS_PER_TASK}"
BUFFER=100000

# Your new flag that triggers: R1 PBAT-like, R2 non-PBAT (after you implement it in Bhmem.java)
NEW_MODE_FLAG="-snm3c"

# RG
RGSM="snm3C"

# -----------------------
# Pick sample for this array task
# -----------------------
prefix="$(awk -v num="${SLURM_ARRAY_TASK_ID}" 'NR==num{print; exit}' "${ACC_LIST}")"
if [[ -z "${prefix}" ]]; then
  echo "ERROR: Could not find prefix for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${ACC_LIST}" >&2
  exit 1
fi

R1="${TRIMDIR}/${prefix}-R1.trimmed.fq.gz"
R2="${TRIMDIR}/${prefix}-R2.trimmed.fq.gz"

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
  echo "ERROR: Missing trimmed FASTQ: ${R1} or ${R2}" >&2
  exit 1
fi

OUTBAM="${OUTDIR}/${prefix}.bhmem.bam"

if [[ -e "${OUTBAM}" ]]; then
  echo "ERROR: Output exists (won't overwrite): ${OUTBAM}" >&2
  exit 1
fi

echo "Sample: ${prefix}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Ref: ${REF_FA}"
echo "Out: ${OUTBAM}"
echo "Threads: ${THREADS}"
echo "Mode flag: ${NEW_MODE_FLAG}"

# -----------------------
# Run Bhmem
# -----------------------
java -Xmx60G \
  -Djava.library.path="${JNI_LIB}" \
  -cp "${CP}" \
  "${MAIN_CLASS}" \
  "${REF_FA}" \
  "${OUTBAM}" \
  "${R1}" \
  "${R2}" \
  -rgId "${prefix}" \
  -rgSm "${RGSM}" \
  -nonDirectional \
  -pbat \
  -buffer "${BUFFER}" \
  -enzymeList "${ENZYME_BED}" \
  -outputMateDiffChr

echo "Bhmem finished: ${OUTBAM}"

# Quick index + basic sanity check
samtools index -@ "${THREADS}" "${OUTBAM}"
samtools flagstat -@ "${THREADS}" "${OUTBAM}" > "${OUTBAM}.flagstat.txt"
echo "Index + flagstat done: ${OUTBAM}.flagstat.txt"
