# snmCAT — Submitted Jobs

| Date | Job name | Job ID | Submit command | Array | Logs | Status |
|---|---|---|---|---|---|---|
| 2026-06-10 | snmCAT_mct | 4307786 | `sbatch codes/02.run_snakemake.sh` | 1-64 | `codes/logs/02.mapping/snakemake.<task>.{out,err}` | submitted (PD) |
| (earlier) | dl_fastq | (see codes/logs/00.download) | `sbatch codes/00.download.sh` | 1-100 | `codes/logs/00.download/dl.<task>.{txt,err}` | done (100 cells) |

## Latest job: 4307786 (yap mct mapping)
- One array task per snakemake Group (64 groups, 100 cells total).
- Each task: `snakemake -d mapping/Group<k> --snakefile .../Snakefile -j 10 ...`
- Outputs per cell under `mapping/Group<k>/`: `allc/<cell>.allc.tsv.gz` (methylation),
  `rna_bam/...feature_count.tsv` (RNA gene counts), `bam/...` (bismark BAMs), `*.stats`.
- After all tasks succeed: run `yap summary` in `mapping/` to build `stats/MappingSummary.csv.gz`.

### Check status
```
squeue -u jmj7858 -j 4307786 -t all
```
### If a task fails
Inspect `codes/logs/02.mapping/snakemake.<task>.err`, fix, then resubmit just that task:
`sbatch --array=<task> codes/02.run_snakemake.sh`
