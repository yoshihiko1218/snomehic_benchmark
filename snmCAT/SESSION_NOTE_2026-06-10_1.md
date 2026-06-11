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

## FINAL OUTCOME (2026-06-10, end of session) — ✅ COMPLETE
Pipeline finished successfully: **100/100 cells mapped**, 0 errors.
- Three bugs found & fixed during monitoring:
  1. **STAR version mismatch** — index `star_2.7.11b` unreadable by env STAR 2.7.3a (`genomeType` FATAL).
     Rebuilt matching index `star_2.7.10a_gencode.v36_sjdb100` (job 4312628; actual builder 2.7.3a).
  2. **Index build crashed in 0s** — `set -euo pipefail` ran before `source ~/.bashrc` (BASHRCSOURCED
     unbound var). Fixed by sourcing bashrc/activating conda BEFORE `set -u`.
  3. **Stale-lock restart** — runner `--unlock` aborted on IncompleteFilesException before clearing
     locks (only 7/64 groups ran). Fixed: `${CMD} --rerun-incomplete --unlock` then `${CMD} --rerun-incomplete`.
- Final job: **4316481** (array 1-64) — all 100 cells, allc + RNA feature_count produced.
- `yap summary -o mapping` → `mapping/stats/MappingSummary.csv.gz` (100×92) + `AllcPaths.tsv` (100).
- QC medians: R1 input 2.52M, R1 map 64.8% / R2 32.6%, mCG 72.8%, mCHH 9.2%, FinalDNA 733k; RNA cols present.
- This MappingSummary.csv.gz is the benchmark QC input (same format as other yap methods per PROJECT_CONTEXT).

## NOMe / HCG-GCH investigation (2026-06-10) — VERDICT: snmCT-seq, NOT NOMe
User asked to extract per-cell HCG & GCH loci, and questioned whether this is snmCAT-seq/snmC2T-seq (NOMe).
- HCG = H-CG (CpG, not preceded by G); GCH = G-CH (GpC accessibility, only real if NOMe/GpC-MTase used). GCG excluded.
- Original ALLC used num_upstr_bases=0 → can't classify HCG/GCH. Recomputed with num_upstr_bases=1
  via `allcools bam-to-allc` from retained dna_reads.bam (codes/05, 06) → mapping/stats/hcg_gch_site_counts.tsv.
- **Code vs paper check** (cemba_data/mapping/mct/mct_bismark_bam_filter.py): yap WITHOUT --nome == paper
  step 4a exactly (XM tag, mCH<=0.5, cov>=3). `--nome` additionally EXCLUDES GpC from the read-level mCH
  (so accessibility reads aren't misbinned as RNA); paper step 5a needs num_upstr_bases=1. yap leaves
  bismark_reference/star_reference/nome_flag_str BLANK in generated files (injection bugs) — all patched.
- **Full --nome re-run** of all 100 cells (job 4329802 → mapping_nome/, config codes/mapping_config_nome.ini,
  runner codes/08, collector codes/09 → mapping_nome/stats/hcg_gch_nome.tsv).
- **RESULT (medians): GCH ≈ HCH background in BOTH pipelines; --nome changed nothing.**
  non-nome: HCG 60.2% / GCH 1.57% / HCH 1.53% ; --nome: HCG 65.4% / GCH 1.53% / HCH 1.59%.
- **VERDICT: this dataset is snmCT-seq (mC + transcriptome), NO NOMe/GpC accessibility.** If NOMe had been
  applied GCH would be 15-40%. The 'mCT' SRA labels are correct; folder name 'snmCAT' is a misnomer.
  HCG (~60-65%) is the real CpG-methylation readout. GCH carries no accessibility signal (only coverage breadth).
- The yap parameter was NOT the cause — proper --nome processing gave identical GCH because there is no GpC methylation.
- Canonical output going forward: mapping_nome/ (cleaner num_upstr_bases=1 ALLC, NOMe-ready context labels).

## FINAL (2026-06-11): found + verified the REAL NOMe data = brain snmC2T-seq
- H1/HEK "mCT" batches (171009 scmCT, 180615 snmCT) both NON-NOMe (GCH ~1.3-1.5%) — confirmed by GEO protocol field + data.
- The NOMe (snmC2T-seq) data in GSE140493 is the BRAIN cortex samples (UMB5577/UMB5580, 190305/190321 mCTseq).
- Downloaded 100 cells (UMB5580, 190321), mapped with full yap mct --nome (mapping_brain/), 0 errors.
- **100-cell medians: HCG 81.2% (CpG meth) / GCH 15.1% (accessibility) / HCH 4.6% (bg); GCH/HCH 3.3x; 92/100 cells GCH>2x HCH.**
- => Real NOMe chromatin accessibility CONFIRMED. Pipeline validated end-to-end. Recipe: RUNBOOK_nome.md (codes 15-19).
