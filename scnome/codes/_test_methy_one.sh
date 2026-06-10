set -uo pipefail
source /home/jmj7858/.bashrc
conda activate scnomehic
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome
BIS=/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/
METHY=05.methy
bam=04.alignment/SRR3729670_1.rmdup.bam
echo "[$(date)] step1 extractor start"
bismark_methylation_extractor -s --ignore 6 --comprehensive --multicore 4 \
    --bedGraph --CX --genome_folder "$BIS" -o "$METHY" "$bam"
echo "[$(date)] step1 done; cov:"
ls -la $METHY/SRR3729670_1.rmdup.bismark.cov.gz
echo "[$(date)] step2 coverage2cytosine --nome-seq start"
coverage2cytosine --nome-seq --CX --genome_folder "$BIS" --dir "$METHY" --gzip \
    -o SRR3729670_1 "$METHY/SRR3729670_1.rmdup.bismark.cov.gz"
echo "[$(date)] step2 done. NOMe outputs:"
ls -la $METHY/SRR3729670_1*NOMe* 2>/dev/null
echo "TESTDONE"
