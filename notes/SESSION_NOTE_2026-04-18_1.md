# Session Note — 2026-04-18 (1)

## Goal
1. Walk-through summary of benchmark folder (done).
2. Check how `summary/trinuc/` is generated and find methylation methods missing.
3. Run bhmem pipeline on scnomehic to obtain per-cell `.trinuc_methy.chrX.txt`.
4. Drop `methylhic_new` from benchmark comparison.

## Findings
- `summary/extract_trinuc.py` + `run_extract_trinuc.sh` pull ACT/ACG/GCT rows from
  per-cell `*.trinuc_methy.chrXX.txt` files. Those files are produced by
  `bissnp_trinuc_sample.pl` (Bis-tools Bis-QC), invoked during `03.bamprocess.sh`
  when Bis-QC runs with `--pattern WCH`.
- Methylation methods previously in `summary/trinuc/`: scnome, smallwood, snmCseq3,
  snmCseq2. Missing: methylhic (files existed, just not extracted), methylhic_new
  (YAP-only, no trinuc files), scnomehic (YAP-only, no trinuc files).
- scnomehic was NOT previously run with bhmem — only YAP (`alignment/` = bowtie2,
  `alignment_bowtie1/`). No `04.bhmem_bam/` existed.

## Actions
- Created `scnomehic/acc_list.txt` — 188 cell prefixes from `scnomehic/fastq/*-R1.fq.gz`.
- Wrote bhmem pipeline scripts in `scnomehic/codes/`:
  - `01.trim.sh`: trim_galore paired, clip 15/5, -j 32, 1-188 array
  - `02.alignment.sh`: Bhmem -nonDirectional -pbat -buffer 100000, 1-188 array
  - `03.bamprocess.sh`: fixmate/markdup/calmd + Bis-QC WCH + Bis-SNP (copied from
    snmCseq3 with paths adjusted), 1-188 array
- Submitted chained SLURM jobs:
  - 5703168 (trim) → 5703169 (bhmem, afterok) → 5703170 (bamprocess, afterok)
- Recorded in `scnomehic/JOBS.md` and `scnomehic/codes/FILES.md`.
- Ran `methylhic` trinuc extraction immediately (59 samples): produced
  `summary/trinuc/methylhic.chr19.txt`.
- Updated `summary/run_extract_trinuc.sh` to add `methylhic` and `scnomehic` lines.
- Removed `methylhic_new` row from `summary/datasets_all.csv`.

## Next steps (after SLURM jobs finish)
1. Verify `scnomehic/04.bhmem_bam/*.calmd.trinuc_methy.chr19.txt` (188 files expected).
2. Run scnomehic line of `run_extract_trinuc.sh` → `summary/trinuc/scnomehic.chr19.txt`.
3. Update `summary/compare_qc_all.py` `load_trinuc_data()` to include
   `methylhic.chr19.txt` and `scnomehic.chr19.txt` (merging NonCG_rate/GCH_rate).
4. Re-run `python compare_qc_all.py --config datasets_all.csv --output all_methods_qc`.

## Files modified
- `scnomehic/acc_list.txt` (new)
- `scnomehic/codes/01.trim.sh` (new)
- `scnomehic/codes/02.alignment.sh` (new)
- `scnomehic/codes/03.bamprocess.sh` (new)
- `scnomehic/codes/FILES.md` (new)
- `scnomehic/JOBS.md` (new)
- `summary/trinuc/methylhic.chr19.txt` (new — 60 rows incl. header)
- `summary/run_extract_trinuc.sh` (added methylhic + scnomehic lines)
- `summary/datasets_all.csv` (dropped methylhic_new row)

## Latest submitted job (per global CLAUDE.md rule)
`5703170` (scnh_bamproc, pending dependency on 5703169). Check with `squeue -u jmj7858`.
