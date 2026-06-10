# snmCAT — Submitted Jobs

| Date | Job name | Job ID | Submit command | Array | Logs | Status |
|---|---|---|---|---|---|---|
| 2026-06-10 | snmCAT_mct (resubmit) | 4314070 | `sbatch codes/02.run_snakemake.sh` | 1-64 | `codes/logs/02.mapping/snakemake.<task>.{out,err}` | submitted (PD) — new 2.7.3a index |
| 2026-06-10 | star_index_2710a | 4312628 | `sbatch codes/04.build_star_index_2.7.10a.sh` | - | `codes/logs/04.star_index/build.{out,err}` | COMPLETED (23 min; builder STAR=2.7.3a, see note) |
| 2026-06-10 | star_index_2710a | 4310364 | `sbatch codes/04.build_star_index_2.7.10a.sh` | - | `codes/logs/04.star_index/build.{out,err}` | **FAILED** (set -u + bashrc unbound var, 0s) |
| 2026-06-10 | snmCAT_mct | 4307786 | `sbatch codes/02.run_snakemake.sh` | 1-64 | `codes/logs/02.mapping/snakemake.<task>.{out,err}` | **CANCELLED** (STAR version mismatch) |
| (earlier) | dl_fastq | (see codes/logs/00.download) | `sbatch codes/00.download.sh` | 1-100 | `codes/logs/00.download/dl.<task>.{txt,err}` | done (100 cells) |

## Issue found & fix (2026-06-10)
- Job 4307786 failed in `rule star`: `FATAL INPUT ERROR: unrecognized parameter name "genomeType"`.
- Cause: index `star_2.7.11b_gencode.v36_sjdb150` built with STAR **2.7.11b**, but yap's `mapping`
  env has STAR **2.7.10a** (cannot read the newer index). Bismark/methylation half was fine.
- NOTE: inside the SLURM job (with `export PATH="$CONDA_PREFIX/bin:$PATH"`), the env's STAR is
  actually **2.7.3a** (conda `mapping/bin/STAR`), not 2.7.10a. The mapping job uses the same
  preamble, so index builder and aligner match (2.7.3a). Dir name "2.7.10a" is cosmetic only.
- Fix: build a matching index `star_2.7.10a_gencode.v36_sjdb100` (job 4310364), repoint
  `mapping_config.ini` + Snakefiles (`codes/03.patch_star_reference.sh` now overwrites the value),
  then resubmit `codes/02.run_snakemake.sh`.

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
