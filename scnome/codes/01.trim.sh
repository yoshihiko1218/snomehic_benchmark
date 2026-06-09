#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-12
#SBATCH --job-name=qc
#SBATCH --output=logs/01.qc_trim/qc_trim.%a.txt
#SBATCH --error=logs/01.qc_trim/qc_trim.%a.err

source /home/jmj7858/.bashrc
cd /gpfs/projects/b1198/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome

conda activate scnomehic

infolder=01.fastq
outfolder=02.fastqc_out
outfolder_fq=03.trimmed_fastq

prefix=`cat acc_list.txt | awk -v num=${SLURM_ARRAY_TASK_ID} 'NR == num'`
echo ${prefix}

# Classify sample by SRR accession number and set trimming accordingly:
#   SRR3729642-SRR3729653 : GM12878 cells -> clip 6 bp from BOTH ends
#                           (paper: "In the case of GM12878 cells 6 bp were
#                            clipped from either end of the read")
#   SRR3729654-SRR3729660 : control / spike-in samples -> excluded from analysis
#   SRR3729661-SRR3729682 : K562 cells    -> clip 6 bp from the 5' end only
srr_num=${prefix#SRR}
if [ "${srr_num}" -ge 3729642 ] && [ "${srr_num}" -le 3729653 ]; then
  celltype=GM12878
  clip_args="--clip_R1 6 --three_prime_clip_R1 6"
elif [ "${srr_num}" -ge 3729654 ] && [ "${srr_num}" -le 3729660 ]; then
  celltype=control
  echo "${prefix} is a control sample (excluded from analysis). Skipping."
  exit 0
else
  celltype=K562
  clip_args="--clip_R1 6"
fi
echo "Cell type: ${celltype}; trim_galore clip args: ${clip_args}"

fastqc --outdir ${outfolder} -t 8 ${infolder}/${prefix}_1.fastq ${infolder}/${prefix}_2.fastq
echo "First fastqc Done"

trim_galore --quality 30 --phred33 --illumina --stringency 1 -e 0.1 ${clip_args} --gzip --length 20 -j 8 -o ${outfolder_fq} --fastqc --fastqc_args "--outdir ${outfolder} -t 8" ${infolder}/${prefix}_1.fastq ${infolder}/${prefix}_2.fastq
echo "Trimming Done"
