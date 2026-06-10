# Session Note 2026-06-10_1 — snmCAT yap pipeline setup

## Goal
Run a yap-based mapping pipeline on the snmCAT dataset, mirroring snmCseq3's
yap workflow (in the parent `benchmark/` folder).

## Investigation findings
- **Dataset**: `snmCAT/fastq/` = 100 single cells, paired-end, SRR10470168.. (GSM4162955..),
  titles `171009_mCT_hs_h1_*` → **snmCT-seq** (methylC + Transcriptome), **human H1 ESC**.
  ~2.2M read pairs/cell, 150bp, already cell-level (demux done by SRA; no inline barcode).
  Download complete (200 fastq.gz = 100 R1 + 100 R2). Meta: `snmCAT/meta.tsv` (8494 SRRs total, 100 chosen).
- **yap**: v1.6.9 at conda env `mapping`. Modes: {mct, mc, m3c, 4m}, optional `--nome`.
- **snmCseq3 yap workflow** (the reference):
  1. `00.rename.sh` — fastq → `{prefix}-R[12].fq.gz`
  2. `yap start-from-cell-fastq -o alignment/ -config codes/mapping_config.ini -fq "fastq/*-R[12].fq.gz"`
     → writes `alignment/snakemake/snakemake_cmd.txt` (one snakemake command per cell-group)
  3. `01.excute_snakemake.sh` — SLURM array, runs Nth line of snakemake_cmd.txt
  - config: snmCseq3 used `mode=m3c`, mm10 refs, barcode V2.
- **hg38 references (all present & verified)**:
  - bismark_ref: `/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark` (has Bisulfite_Genome)
  - genome_fasta: `.../hg38_bismark/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa` (or hg38/ same file)
  - star_ref: `.../hg38/star_2.7.11b_gencode.v36_sjdb150` (SA, SAindex, Genome present)
  - gtf: `.../hg38/gencode.v36.annotation.gtf`
  - chrom_size: `.../hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes`

## Decision
- Mode: **mct** (matches snmCT data; separates DNA/RNA before mC calling). snmCT is NOT NOMe → no `--nome`.
- Awaiting user confirmation on mode + cell count before executing.

## Actions taken
1. `codes/01.rename_symlink.sh` — created 100×2 yap symlinks `SRR*-R{1,2}.fq.gz`. Ran locally.
2. Generated `codes/mapping_config.ini`:
   `yap default-mapping-config --mode mct -v V2 --bismark_ref <hg38_bismark> --genome_fasta <...fa> --star_ref <hg38/star_2.7.11b_gencode.v36_sjdb150> --gtf <gencode.v36.annotation.gtf> --chrom_size_path <...chrom.sizes>`
3. `yap start-from-cell-fastq -o mapping -config codes/mapping_config.ini -fq "fastq/*-R[12].fq.gz"`
   → 64 Group Snakefiles + `mapping/snakemake/snakemake_cmd.txt` (64 commands, 100 cells). Confirmed bismark+STAR engine.
4. **Found yap 1.6.9 bug**: Snakefile uses `{bismark_reference}` and `{star_reference}` but yap doesn't inject them
   → snakemake `NameError`. Fixed via `codes/03.patch_star_reference.sh` (injects both into all 64 Snakefiles). Idempotent.
5. Validated with `snakemake -n` (dry-run) on Group45 → full 18-rule DAG builds cleanly, no errors.
6. Submitted `codes/02.run_snakemake.sh` → **job 4307786** (array 1-64), PD. Logs `codes/logs/02.mapping/`.

## Next steps
- Monitor 4307786 (squeue). On failure: read `codes/logs/02.mapping/snakemake.<task>.err`, fix, resubmit that task.
- After all succeed: `yap summary` in `mapping/` → `stats/MappingSummary.csv.gz` for the benchmark QC.
- `mapping_probe1/` is a throwaway probe dir (can be deleted).

## Decisions confirmed by user
- Mode = **mct** (methylation + RNA). Scope = **all 100 cells**.
