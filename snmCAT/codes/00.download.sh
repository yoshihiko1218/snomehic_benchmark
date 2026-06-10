#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100%10
#SBATCH --job-name=dl_fastq
#SBATCH --output=logs/00.download/dl.%a.txt
#SBATCH --error=logs/00.download/dl.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT

outfolder=fastq
listfile=codes/download_list.txt

# this task's cell: semicolon-separated R1;R2 ENA ftp paths
line=`cat ${listfile} | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo "Task ${SLURM_ARRAY_TASK_ID}: ${line}"

# download each read with resume (-c) + retries
echo "${line}" | tr ';' '\n' | while read url; do
    [ -z "${url}" ] && continue
    echo "Downloading https://${url}"
    wget -c -t 5 --waitretry=10 --read-timeout=60 -P ${outfolder} "https://${url}"
done
echo "Download Done for task ${SLURM_ARRAY_TASK_ID}"
