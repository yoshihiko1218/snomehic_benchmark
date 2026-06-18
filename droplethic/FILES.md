# FILES.md — droplethic/

Droplet Hi-C dataset. Processed via a nested pipeline project rather than the
numbered-stage convention used by other datasets.

| Entry | What it is |
|-------|-----------|
| `acc_list.txt` | SRA accession list. |
| `fastq/` | Raw FASTQ (gitignored). |
| `my_project/` | Pipeline working dir. Per-cell QC table read by `summary/qc.ipynb`: `my_project/SRR27586278_hg38.per_cell_qc.valid.tsv`. |

Note: `summary/qc.ipynb` loads droplet Hi-C QC from `droplethic/my_project/*.per_cell_qc.valid.tsv`.
