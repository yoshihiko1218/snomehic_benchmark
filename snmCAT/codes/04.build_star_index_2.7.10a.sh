#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 08:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=star_index_2710a
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/04.star_index/build.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/04.star_index/build.err

# Build a STAR index with the SAME STAR version yap uses (2.7.10a in the `mapping`
# env). The existing star_2.7.11b_* index is unreadable by 2.7.10a
# (FATAL: unrecognized parameter name "genomeType"). hg38 + gencode.v36, sjdbOverhang 100.

# NOTE: source bashrc / activate conda BEFORE `set -u`; ~/.bashrc references an
# unbound var (BASHRCSOURCED) that would abort the job under set -u.
source /home/jmj7858/.bashrc
conda activate mapping
export PATH="$CONDA_PREFIX/bin:$PATH"

set -euo pipefail

echo "STAR: $(which STAR)"; STAR --version

REF=/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38
FASTA=${REF}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
GTF=${REF}/gencode.v36.annotation.gtf
OUT=${REF}/star_2.7.10a_gencode.v36_sjdb100

mkdir -p "${OUT}"
STAR --runMode genomeGenerate \
     --genomeDir "${OUT}" \
     --genomeFastaFiles "${FASTA}" \
     --sjdbGTFfile "${GTF}" \
     --sjdbOverhang 100 \
     --runThreadN "${SLURM_CPUS_PER_TASK:-12}"

echo "[$(date)] STAR index build DONE -> ${OUT}"
ls -la "${OUT}"
