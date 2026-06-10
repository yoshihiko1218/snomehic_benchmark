curl -s 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA589980&result=read_run&fields=run_accession,experiment_alias,experiment_title,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp&format=tsv' \
  | tail -n +2 \
  | sort -t$'\t' -k2,2 \
  | head -100 \
  | cut -f4 > download_list.txt

sbatch 00.download.sh
