# Session Note — 2026-05-28 (1)

## Goal
Rerun the scNOME-HiC processing using the UPDATED Snakemake pipeline at
`/gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline`.

User decisions (via clarifying questions):
- **Workdir:** `gm_sc_new` (hg38 GM12878 cellline) =
  `/gpfs/projects/b1198/epifluidlab/yoshii/sc_nomehic_cellline/gm_sc_new`.
  This is where the pipeline `configs/config.yaml` already points (start_from: raw).
- **Scope:** trimming → alignment → bam processing → methylation. **NO hicluster.**
  (hicluster is already excluded from the pipeline's default `rule all`.)

## Pipeline state observed
- Pipeline git: branch master, recent commits (to 2026-05-27) all on hicluster
  + gzip-all-text-outputs (2026-05-18) + bisqc autosome fallback (2026-05-15)
  + enzyme_list config (2026-05-15). One uncommitted edit: `concatcells`
  mem_mb=64000 (hicluster only — irrelevant to this run).
- Default `rule all` targets (raw mode): 02.fastqc, 04.alignment summary.txt.gz,
  07.bistools qc hist.txt.gz, 07.bistools methylation 6plus2.bed.gz,
  08.methylation *_methylation.txt.gz. No hicluster.
- gm_sc_new already has prior outputs: 04.alignment_snakemake, 08.methylation_snakemake,
  GCH/HCG merged, etc. So this is a RESUME, not a from-scratch run.

## Actions
- (in progress) Running `snakemake -np` dry-run to see what is missing/stale.

## Dry-run result (resume scope, before the from-scratch decision)
With `--rerun-triggers mtime`: 1881 jobs would run (qc 188, bistools 188,
methylation 1504). Trim/mapping/bamprocess were already up-to-date.

## User pivot → FROM SCRATCH
User: "redo from scratch because I fixed alignment and bistools step."
So existing outputs are stale → full `--forceall` rerun of the whole DAG
(demux/trim → mapping → bamprocess → qc → bisqc → bistools → methylation),
188 cells. No hicluster (not in default rule all).

Decisions:
- Mode: `--forceall`, overwrite in place (keep top-level merged files).
- Launch: dedicated sbatch DRIVER job (persistent).

## Pre-flight fixes
- DpnII enzyme files were gzipped (May 17 sweep). mapping rule's `-enzymeList`
  fallback `{reference}.DpnII.span_region.bedgraph` was missing → `gunzip -k`
  restored the uncompressed file (kept .gz).
- Confirmed the already-running hicluster driver (PID 3327535) is on workdir
  `sc_nomehic/snakemake`, NOT gm_sc_new → no collision; gm_sc_new lock empty.
- Partition: driver on `genomicslong` (10-day max); `genomics` caps at 48h
  which fits each child rule. b1042 has GCC partition access.

## Files created/modified
- `codes/dryrun_pipeline.sh` (new) — pipeline dry-run helper.
- `codes/run_pipeline_gm_sc_new.sh` (new) — driver sbatch job.
- `codes/FILES.md`, `JOBS.md` — updated with the rerun + driver job.
- reference: gunzipped `.DpnII.span_region.bedgraph`.

## Submitted job (per global CLAUDE.md rule)
**9439299** `scnh_driver` (genomicslong, 7d). PENDING (Priority) at submit.
Driver log: `logs/pipeline_rerun/driver.9439299.{out,err}`.
Child jobs appear under `squeue -u jmj7858` (UUID names, partition genomics).

## Next
- When driver starts, confirm via driver.9439299.out that `sbatch` works on the
  node and snakemake begins submitting child jobs (mapping first after trim).
- Watch for mapping failures (enzyme_list path) and bistools failures.
- Check progress: `squeue -u jmj7858`.
