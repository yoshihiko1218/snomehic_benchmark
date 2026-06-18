# FILES.md — methylhic/

Methyl-HiC dataset (numbered-stage pipeline).

| Entry | What it is |
|-------|-----------|
| `acc_list.txt` | SRA accession list. |
| `codes/` | Pipeline scripts (trim, align, methylation, QC). |
| `01.fastq/` | Raw FASTQ (gitignored). |
| `02.fastqc_out/` | FastQC reports (gitignored). |
| `03.trimmed_fastq/` | Trimmed FASTQ (gitignored). |
| `04.alignment/` | Aligned BAMs + summaries (BAMs/large txt gitignored). |
| `alignment/` | Snakemake per-Group alignment working dir (`Group*/`, gitignored data). |
| `fastq_renamed/` | Renamed FASTQ inputs for the pipeline. |
| `logs/` | Job logs. |

Note: two alignment dirs exist (`04.alignment/` and `alignment/Group*/`) from
different pipeline stages; both retained because Snakefiles reference these paths.
