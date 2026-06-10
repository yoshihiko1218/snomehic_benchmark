# FILES — scnome/codes

Pipeline for the scNOMe-seq (Pott 2017) benchmark dataset (Trim Galore + Bismark SE,
each mate processed independently, hg38). Scripts are named by execution stage:
`01` trim → `02`/`02k` align-or-merge + dedup → `03`/`03k` methylation → `04` QC.
A `k` suffix marks the K562 variant; `submit_*` are chain drivers; `util_*` are
one-off maintenance helpers; `*.py` are helpers called by the stage scripts.

## Cell-type roster (`../acc_list.txt`, 34 cells after control removal)
- **GM12878** — SRR3729642–SRR3729653 → lines 1–12. **1 run per cell.**
- **K562** — SRR3729661–SRR3729682 → lines 13–34. **2 runs per cell** (consecutive
  SRR pairs); 11 merged cells `K562_01`..`K562_11` in `../acc_list_k562_cells.tsv`
  (ids in `../acc_list_k562_cellids.txt`).
- **Controls** (spike-in, EXCLUDED) — SRR3729654–SRR3729660 → removed from
  `acc_list.txt`; original list preserved as `../acc_list_with_controls.txt`.

## Pipeline stages (SLURM array jobs)
- **01.trim.sh** — FastQC + Trim Galore over `acc_list.txt` runs. Per-SRR clip:
  GM12878 → `--clip_R1 6 --three_prime_clip_R1 6` (6 bp BOTH ends, Pott 2017);
  K562 → `--clip_R1 6` (5' only); controls → skipped. Inputs are `01.fastq/*.fq.gz`.
- **02.align_dedup.sh** — *(GM)* Bismark align (`--non_directional --score_min L,0,-0.2`,
  hg38) each mate → `samtools sort | markdup` → `.rmdup.bam` → `addreplacerg`
  → `.rmdup.RG.bam` (+index) → `bam_summary_universal.py`. markdup WITHOUT `-r`.
- **02k.merge_dedup_k562.sh** — *(K562)* per cell: `samtools cat` the two runs'
  per-mate bismark BAMs → sort → markdup (removes within+cross-run dups) →
  `K562_NN_<m>.rmdup(.RG).bam` → summary. Reuses existing per-run bismark BAMs.
- **03.methy_extract.sh** — *(GM)* Bismark NOMe extraction: `bismark_methylation_extractor
  -s --ignore 6 --bedGraph --CX` then `coverage2cytosine --nome-seq` → per-mate
  `.NOMe.{CpG,GpC}.cov.gz`. Resume-safe.
- **03k.methy_extract_k562.sh** — *(K562)* same as 03 but over the merged cell list.

## Stage 04 — QC (each is an array job over a LISTFILE)
- **04.qc_per_cell.sh** — per-cell QC via `scnome_qc_per_cell.py` → `qc_stats/<cell>.qc_stats.csv`
  (trim, Bismark, BAM, methyl, site & trinuc columns).
- **04.qc_nome_sites.sh** — `nome_qc_sites_trinuc.py`: detected HCG/GCH site counts
  (from NOMe cov) + chrM/chr21 trinuc proxy → `qc_stats/<cell>.nome_qc.tsv`.
- **04.qc_bissnp_trinuc.sh** — ORIGINAL BisSNP trinuc QC (BisulfiteGenotyper,
  `-minPatConv 0.8`, dbSNP-aware) on chrM+chr21 → `04.alignment/<cell>_<m>.rmdup.RG.trinuc_methy.{chrM,chr21}.txt`.
  **Requires Java 8** (`module load java/jdk1.8.0_191`; GATK-3.8 fails on Java 21).
- **04.qc_collect.sh** — per-cell QC over `acc_list.txt` then aggregate via
  `collect_scnome_qc.py` → `../scnome_qc_summary.csv`.

## Chain drivers
- **submit_gm.sh** — GM chain: 01.trim → 02.align_dedup → 03.methy_extract → 04.qc_per_cell (afterok).
- **submit_k562.sh** — K562 chain: 02k.merge_dedup_k562 → 03k.methy_extract_k562 (afterok).

## Helpers (`*.py`, called by stages)
- **scnome_qc_per_cell.py** — per-cell QC metrics parser.
- **collect_scnome_qc.py** — aggregate per-cell QC CSVs into one summary.
- **nome_qc_sites_trinuc.py** — Bismark-native HCG/GCH site counts + trinuc proxy.
- **se_bam_summary.py** — single-end BAM mapping summary (alternate; not in core run).

## Utilities (one-off, DRY-RUN by default; `--yes` to act)
- **util_clear_gm_stale.sh** — delete STALE GM12878 outputs (old 5'-only trimming).
- **util_clear_control_files.sh** — delete all control-sample files (SRR3729654–3729660).

## Tracking
- **JOBS.md** — submitted-job tracking. **qc.ipynb** — QC exploration notebook.
- **_test_c2c.sh / _test_methy_one.sh** — ad-hoc test scripts (not pipeline steps).
