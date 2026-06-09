# FILES — scnome/codes

Pipeline for the scNOMe-seq (Pott 2017) benchmark dataset (Trim Galore + Bismark SE,
each mate processed independently, hg38). Sample roster is `../acc_list.txt`.

## Cell-type roster (`../acc_list.txt`, 34 cells after control removal)
- **GM12878** — SRR3729642–SRR3729653  → lines 1–12
- **K562**    — SRR3729661–SRR3729682  → lines 13–34
- **Controls** (spike-in, EXCLUDED)  — SRR3729654–SRR3729660  → removed from
  `acc_list.txt`; original list preserved as `../acc_list_with_controls.txt`.

## Pipeline scripts (SLURM array jobs; arrays currently scoped to GM = 1-12)
- **01.trim.sh** — FastQC + Trim Galore. Classifies each sample by SRR number and
  clips accordingly: GM12878 → `--clip_R1 6 --three_prime_clip_R1 6` (6 bp BOTH
  ends, per Pott 2017); K562 → `--clip_R1 6` (5' only); controls → skipped.
- **02.alignment.sh** — Bismark align (`--non_directional --score_min L,0,-0.2`,
  hg38) of each mate → `samtools sort | markdup` → `.rmdup.bam` → `addreplacerg`
  → `.rmdup.RG.bam` (+index) → `bam_summary_universal.py`. markdup used without
  `-r` (identical dedup to all cell types; only trimming differs for GM).
- **03.methy_extract.sh** — Bismark NOMe extraction: `bismark_methylation_extractor
  -s --ignore 6 --bedGraph --CX` then `coverage2cytosine --nome-seq` → per-mate
  `.NOMe.{CpG,GpC}.cov.gz`. Resume-safe.
- **run_qc.sh** — array job: per-cell QC via `scnome_qc_per_cell.py` → `qc_stats/`.
- **run_qc_and_collect.sh** — single job: per-cell QC over all of `acc_list.txt`
  then aggregate via `collect_scnome_qc.py` → `../scnome_qc_summary.csv`.

## Helper / analysis scripts
- **scnome_qc_per_cell.py** — per-cell QC metrics parser.
- **collect_scnome_qc.py** — aggregate per-cell QC CSVs into one summary.
- **se_bam_summary.py** — single-end BAM mapping summary (alternate; not in core run).
- **clear_gm_stale.sh** — DRY-RUN by default; `--yes` deletes STALE GM12878 outputs
  (SRR3729642–3729653) built with the old 5'-only trimming, so the rerun
  regenerates them. K562/controls untouched.
- **clear_control_files.sh** — DRY-RUN by default; `--yes` deletes all control
  sample files (SRR3729654–3729660) across data dirs.
- **qc.ipynb** — QC exploration notebook.
- **JOBS.md** — submitted-job tracking.

## GM12878 rerun order (after `--three_prime_clip_R1 6` fix)
1. `bash codes/clear_gm_stale.sh --yes`   (clear stale GM outputs)
2. `sbatch codes/01.trim.sh`              (array 1-12)
3. `sbatch --dependency=afterok:<j1> codes/02.alignment.sh`
4. `sbatch --dependency=afterok:<j2> codes/03.methy_extract.sh`
5. `sbatch --dependency=afterok:<j3> codes/run_qc.sh`  (then run_qc_and_collect / collect)
