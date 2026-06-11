#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100%10
#SBATCH --job-name=dl_brain
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/15.download_brain/dl.%a.txt
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/15.download_brain/dl.%a.err

# Download 100 cells from 190321_mCTseq_hs_29yr (UMB5580) = snmC2T-seq (NOMe) batch.
source /home/jmj7858/.bashrc
BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
outfolder="${BASE}/fastq_brain"
listfile="${BASE}/codes/download_list_brain.txt"
mkdir -p "${outfolder}" "${BASE}/codes/logs/15.download_brain"

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${listfile}")
echo "Task ${SLURM_ARRAY_TASK_ID}: ${line}"
echo "${line}" | tr ';' '\n' | while read url; do
    [ -z "${url}" ] && continue
    echo "Downloading https://${url}"
    wget -c -t 5 --waitretry=10 --read-timeout=60 -P "${outfolder}" "https://${url}"
done
echo "Download Done for task ${SLURM_ARRAY_TASK_ID}"
