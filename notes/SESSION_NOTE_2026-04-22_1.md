# Session Note — 2026-04-22 (1)

## Context
Continuing scnomehic trinuc effort. Second in-benchmark bhmem run (hg38) also
produced 0 mapped reads. Diagnosis: read-name suffix `_1`/`_2` on R1/R2 breaks
BWA pairing.

## User guidance
User said the files in `scnomehic/fastq/` are renamed; pointed at
`/gpfs/projects/b1198/epifluidlab/yoshii/meningioma_scnomehic/data_batch_1` for
pipeline reference, noting only reference paths have changed.

## Findings
1. `meningioma_scnomehic/data_batch_1/00-07*.sh` is the canonical scNOME-HiC
   pipeline: demultiplex → fastqc/trim (trim_galore, clip 15/5) → bhmem
   (`-nonDirectional -pbat -buffer 1000`) → bamprocess → Bis-QC.
2. Original raw scnomehic fastqs with proper Illumina headers live at
   `/projects/b1198/epifluidlab/yoshii/sc_nomehic_cellline/gm_sc_new/00.raw_data/`
   (47 libraries, `{prefix}_R{1,2}_001.fastq.gz`).
3. **`sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake/` already contains full
   bhmem output for all 188 cells** — `.calmd.bam` (~19M mapped reads),
   `.calmd.trinuc_methy.chr21.txt`, `.calmd.trinuc_methy.chrM.txt`. Same pipeline
   as meningioma, just different reference paths (per user note).
4. Therefore no bhmem re-run is needed in the benchmark folder — just extract.

## Actions
- Updated `summary/run_extract_trinuc.sh`: scnomehic line now points at
  `sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake/`.
- Ran `extract_trinuc.py` → `summary/trinuc/scnomehic.chr21.txt` (188 cells).
  Values look healthy: noncpg ~1.1% (good BS conversion), endo ~40% (CpG
  methylation, typical GM12878), exo ~65–73% (NOMe GCH accessibility).
- Updated `PROJECT_CONTEXT.md`: scNOME-HiC bhmem source now points to
  sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake.
- Rewrote `scnomehic/JOBS.md` with full failure post-mortem and cleanup list.

## Files modified today
- `summary/run_extract_trinuc.sh` (scnomehic folder path)
- `summary/trinuc/scnomehic.chr21.txt` (new, 188 cells)
- `PROJECT_CONTEXT.md` (scNOME-HiC bhmem path note)
- `scnomehic/JOBS.md` (full post-mortem)

## Next steps
1. Update `summary/compare_qc_all.py` `load_trinuc_data()` to include
   `methylhic.chr19.txt` and `scnomehic.chr21.txt` (chr21 — not chr19!).
2. Re-run `compare_qc_all.py --config datasets_all.csv --output all_methods_qc`.
3. Inspect `all_methods_qc.methylation.pdf` to verify NonCG_rate + GCH_rate now
   populated for scNOME-HiC and methylhic (and confirm methylhic_new is gone).
4. Decide whether to delete large in-benchmark waste:
   - `scnomehic/03.trimmed_fastq/` (~150 GB)
   - `scnomehic/04.bhmem_bam/` (hg38 broken)
   - `scnomehic/04.bhmem_bam_mm10_wrong/` (mm10 broken)

## Outstanding pending actions
- Cleanup of broken bhmem bam dirs (awaiting user confirmation).
- Update `compare_qc_all.py` trinuc merging.
