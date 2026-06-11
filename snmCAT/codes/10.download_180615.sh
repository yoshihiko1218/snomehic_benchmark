#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100%10
#SBATCH --job-name=dl_180615
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/10.download_180615/dl.%a.txt
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/10.download_180615/dl.%a.err

# Download 100 cells from the 180615_mCT_hs_h1-hek293 NOMe (snmCAT-seq) batch.
source /home/jmj7858/.bashrc
BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
outfolder="${BASE}/fastq_180615"
listfile="${BASE}/codes/download_list_180615.txt"
mkdir -p "${outfolder}" "${BASE}/codes/logs/10.download_180615"

# this task's cell: semicolon-separated R1;R2 ENA ftp paths (Nth line)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${listfile}")
echo "Task ${SLURM_ARRAY_TASK_ID}: ${line}"

echo "${line}" | tr ';' '\n' | while read url; do
    [ -z "${url}" ] && continue
    echo "Downloading https://${url}"
    wget -c -t 5 --waitretry=10 --read-timeout=60 -P "${outfolder}" "https://${url}"
done
echo "Download Done for task ${SLURM_ARRAY_TASK_ID}"
