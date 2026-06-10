source /home/jmj7858/.bashrc
conda activate scnomehic
set -o pipefail
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome
BIS=/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/
METHY=05.methy
echo "[$(date)] c2c --nome-seq (no --CX) start"
coverage2cytosine --nome-seq --genome_folder "$BIS" --dir "$METHY" --gzip \
    -o SRR3729642_1 "$METHY/SRR3729642_1.rmdup.bismark.cov.gz"
echo "[$(date)] exit=$?"
echo "=== outputs ==="
ls -la $METHY/SRR3729642_1* | grep -iE "nome|gpc|cpg"
echo "TESTDONE"
