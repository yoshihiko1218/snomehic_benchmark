# FILES.md — methylhic_new/

Methyl-HiC newer batch (Snakemake per-Group alignment; flat layout).

| Entry | What it is |
|-------|-----------|
| `acc_list.txt` | SRA accession list. |
| `codes/` | Pipeline scripts. |
| `fastq/` | Raw FASTQ (gitignored). |
| `fastq_renamed/` | Renamed FASTQ inputs for the pipeline. |
| `alignment/` | Snakemake per-Group alignment working dir (`Group*/` with `Snakefile`; bam/fastq/allc/hic gitignored). |
| `logs/` | Job logs. |

Note: `summary/qc.ipynb` compares methylhic_new vs scnomehic; alignment Snakefiles
hardcode this folder's absolute paths, so its internal layout is left unchanged.
