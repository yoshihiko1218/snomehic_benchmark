# snmCAT — Submitted Jobs

## ACTIVE: re-run on CORRECT NOMe batch 180615_mCT_hs_h1-hek293 (snmCAT-seq)
Old 171009/171101 (scmCT-seq, non-NOMe) data DELETED by user. Now processing the GpC-treated batch.
| Step | Job/script | Status |
|---|---|---|
| download 100 cells → fastq_180615/ | 4342733 (`codes/10.download_180615.sh`) | running/queued |
| rename symlinks | `codes/11.rename_symlink_180615.sh` | after download |
| yap start-from-cell-fastq (--nome config) → mapping_180615/ | `codes/mapping_config_nome.ini` | after rename |
| patch Snakefiles (refs + nome_flag) | `codes/13.patch_nome_180615.sh` | after start |
| map (yap mct --nome) | `codes/12.run_180615_nome.sh` (`sbatch --array=1-N`) | after patch |
| summary + HCG/GCH | `yap summary -o mapping_180615`; adapt `codes/09` | after map |
Expected if NOMe is real: GCH_mc_rate ~15-40% >> HCH (vs ~1.5% in the old non-NOMe batch).


| Date | Job name | Job ID | Submit command | Array | Logs | Status |
|---|---|---|---|---|---|---|
| 2026-06-10 | snmCAT_mct (resubmit2) | 4316481 | `sbatch codes/02.run_snakemake.sh` | 1-64 | `codes/logs/02.mapping/snakemake.<task>.{out,err}` | ✅ **COMPLETE** — 100/100 cells, 0 errors |
| 2026-06-10 | snmCAT_mct (resubmit) | 4314070 | `sbatch codes/02.run_snakemake.sh` | 1-64 | (same) | partial: 7 groups OK (allc+RNA confirmed), rest blocked by stale-lock/unlock bug |
| 2026-06-10 | star_index_2710a | 4312628 | `sbatch codes/04.build_star_index_2.7.10a.sh` | - | `codes/logs/04.star_index/build.{out,err}` | COMPLETED (23 min; builder STAR=2.7.3a, see note) |
| 2026-06-10 | star_index_2710a | 4310364 | `sbatch codes/04.build_star_index_2.7.10a.sh` | - | `codes/logs/04.star_index/build.{out,err}` | **FAILED** (set -u + bashrc unbound var, 0s) |
| 2026-06-10 | snmCAT_mct | 4307786 | `sbatch codes/02.run_snakemake.sh` | 1-64 | `codes/logs/02.mapping/snakemake.<task>.{out,err}` | **CANCELLED** (STAR version mismatch) |
| (earlier) | dl_fastq | (see codes/logs/00.download) | `sbatch codes/00.download.sh` | 1-100 | `codes/logs/00.download/dl.<task>.{txt,err}` | done (100 cells) |

## NOMe investigation (2026-06-10) — VERDICT: snmCT-seq (no NOMe)
| Job | What | Status |
|---|---|---|
| 4329802 | full 100-cell `--nome` re-run (`codes/08`, → `mapping_nome/`) | ✅ COMPLETE, 0 errors |
| 4340573 | NOMe HCG/GCH/HCH rate collector (`codes/09`) | running → `mapping_nome/stats/hcg_gch_nome.tsv` |
- Medians: non-nome HCG 60.2% / GCH 1.57% / HCH 1.53% ; --nome HCG 65.4% / GCH 1.53% / HCH 1.59%.
- GCH ≈ HCH background in both; `--nome` changed nothing → **no GpC/NOMe signal → snmCT-seq, not snmCAT-seq.**
- HCG = real CpG methylation. Tables: `mapping/stats/hcg_gch_site_counts.tsv` (non-nome),
  `mapping_nome/stats/hcg_gch_nome.tsv` (nome). Canonical output: `mapping_nome/` (num_upstr_bases=1 ALLC).

## ✅ PIPELINE COMPLETE (2026-06-10)
- All 100 cells mapped (yap mct: bismark methylation + STAR/featureCounts RNA), 0 errors.
- `yap summary -o mapping` built **`mapping/stats/MappingSummary.csv.gz`** (100 cells × 92 cols)
  and `mapping/stats/AllcPaths.tsv` (100 ALLC paths). This is the benchmark QC input.
- Per-cell outputs: `mapping/Group*/allc/<cell>.allc.tsv.gz` (methylation),
  `mapping/Group*/rna_bam/*.feature_count.tsv` (RNA gene counts).
- QC medians: R1InputReads 2.52M, R1 MappingRate 64.8% (R2 32.6%), mCG 72.8%, mCHH 9.2%,
  FinalDNAReads 733k. RNA cols present (GenesDetected, RNA/(DNA+RNA), etc.).
- Three bugs fixed en route: (1) STAR 2.7.11b index vs 2.7.3a aligner → rebuilt index;
  (2) `set -u` before `source ~/.bashrc` → reorder; (3) `--unlock` aborting on
  IncompleteFilesException → add `--rerun-incomplete` to the unlock call.

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
