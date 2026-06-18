# FILES.md — nagano/

Nagano et al. single-cell Hi-C (scHi-C) dataset.

| Entry | What it is |
|-------|-----------|
| `acc_list.txt` | SRA accession list for the cells used. |
| `codes/` | Pipeline scripts (download, trim, align, QC). |
| `fastq/` | Raw FASTQ (gitignored). |
| `trimmed_fastq/` | Adapter/quality-trimmed FASTQ (gitignored). |
| `alignment/` | Aligned BAMs + per-cell `.summary.txt` (BAMs gitignored). |
| `qc_stats/` | Per-cell QC tables. |
| `qc_stats_test/` | Scratch/test QC outputs. |
| `qc_summary.csv` | Aggregated per-cell QC summary. |
| `schic_qc_per_cell.py`, `schic_qc.py`, `collect_qc.py`, `run_qc.sh` | QC computation + aggregation. |
| `fastp.html` / `fastp.json` | fastp trimming reports. |
| `RESULTS.md` | Results notes. |
| `logs/` | Job logs. |
