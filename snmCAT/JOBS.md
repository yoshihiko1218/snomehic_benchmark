# snmCAT — Submitted Jobs

## ACTIVE: re-run on CORRECT NOMe batch 180615_mCT_hs_h1-hek293 (snmCAT-seq)
Old 171009/171101 (scmCT-seq, non-NOMe) data DELETED by user. Now processing the GpC-treated batch.
| Step | Job/script | Status |
|---|---|---|
| download 100 cells → fastq_180615/ | 4342733 | ✅ done (200 files, 0 errors) |
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

## CORRECTION (2026-06-11): 180615 H1/HEK is ALSO non-NOMe
- 180615 batch mapped (job 4345601, --nome, mapping_180615/, 100 cells). GCH median **1.3%** ≈ HCH background → NOT NOMe.
- Earlier ID of 180615 as NOMe was WRONG: its GEO **protocol field = "snmCT-seq"** (no "2"); the GpC-MTase
  text was the generic series extraction blob. 171009=scmCT-seq, 180615=snmCT-seq → both non-NOMe (confirmed by data).
- The ONLY NOMe (snmC2T-seq) data in GSE140493 is the **brain UMB5577/UMB5580** samples (190305/190321 mCTseq,
  GSM4167187 protocol field="snmC2T-seq", ~4372 cells). These are human cortex, NOT H1. Not yet verified with data.

## BRAIN snmC2T-seq (REAL NOMe) — 190321_mCTseq UMB5580
| Step | Job/script | Status |
|---|---|---|
| download 100 cells → fastq_brain/ | 4369950 (codes/15) | ✅ done (200 files, 0 errors) |
| map yap mct --nome → mapping_brain/ | 4372418 (codes/18, array 1-64) | running |
| summary + GCH | yap summary -o mapping_brain; codes/19 | after map |
Decisive test: median GCH_mc_rate — expect ~15-40% if real NOMe (vs 1.5% / 1.3% for the H1 batches).

## ✅ NOMe CONFIRMED on brain snmC2T-seq (2026-06-11)
Early preview (map 4372418, mapping_brain/): GCH ELEVATED — SRR10471412 HCG=85.1%/GCH=19.1%/HCH=3.8%;
SRR10470776 HCG=87.1%/GCH=15.8%/HCH=5.2%. GCH ~16-19% >> HCH ~4-5% = real NOMe accessibility.
vs H1 171009 (GCH 1.5%) and 180615 (GCH 1.3%) = non-NOMe. The --nome pipeline works; brain UMB5580 IS the NOMe data.

## ✅ FINAL: brain snmC2T-seq NOMe complete (100 cells, job 4372418 + collector 4396577)
Medians (mapping_brain/stats/hcg_gch_nome.tsv): **HCG 81.2% / GCH 15.1% / HCH 4.6%** ; GCH/HCH = 3.3x ;
92/100 cells GCH>2x HCH. MappingSummary.csv.gz + AllcPaths.tsv built (yap summary -o mapping_brain).
=> Real NOMe accessibility confirmed. vs non-NOMe H1: 171009 GCH 1.5%, 180615 GCH 1.3%.
mapping_brain/ is the canonical snmCAT/NOMe dataset. See RUNBOOK_nome.md.
