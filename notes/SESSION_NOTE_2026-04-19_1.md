# Session Note — 2026-04-19 (1)

## Context carried over from 2026-04-18
Ran bhmem pipeline on scnomehic (188 cells) against mm10 — all 3 stages completed
but produced 0 mapped reads per cell because **scNOME-HiC data is hg38 (human),
not mm10**. PROJECT_CONTEXT.md had wrongly listed scNOME-HiC as mm10. Confirmed
hg38 via YAP `mapping_config.ini` (bismark_reference = hg38_bismark, reference_fasta
= GRCh38) and via YAP BAM header (chr1=248956422 = hg38).

## Fix applied today
1. Moved failed output: `04.bhmem_bam → 04.bhmem_bam_mm10_wrong` (kept for post-mortem).
2. Updated `codes/02.alignment.sh`:
   - `REF_FA = /projects/b1198/epifluidlab/yoshii/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa`
   - `ENZYME_BED = .../GCA_000001405.15_GRCh38_no_alt_analysis_set.DpnII.span_region.bedgraph`
   - Log dir → `logs/02.bhmem_hg38/`
3. Updated `codes/03.bamprocess.sh`:
   - `reference = .../hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set`
   - `vcf = .../hg38/Homo_sapiens_assembly38.dbsnp138.vcf`
   - Resume check now looks for `chr21` not `chr19` trinuc
   - Log dir → `logs/03.bamprocess_hg38/`
4. Resubmitted stage 2 + 3 (skipped stage 1 since trimmed FASTQs are genome-agnostic):
   - `5776407` (scnh_bhmem, hg38) — array 1-188
   - `5776408` (scnh_bamproc, hg38) — array 1-188, `afterok:5776407`
5. Updated `summary/run_extract_trinuc.sh`: scnomehic suffix `.chr19.txt → .chr21.txt`,
   output renamed accordingly.
6. Updated `PROJECT_CONTEXT.md`: scNOME-HiC now listed as hg38, 188 cells, GM12878;
   also dropped the methylhic_new row.

## Files modified today
- `scnomehic/codes/02.alignment.sh` (ref paths + log dir)
- `scnomehic/codes/03.bamprocess.sh` (ref/vcf paths, chr21 resume check, log dir)
- `scnomehic/04.bhmem_bam → 04.bhmem_bam_mm10_wrong` (renamed)
- `scnomehic/JOBS.md` (updated with failed + re-run records)
- `summary/run_extract_trinuc.sh` (scnomehic → chr21)
- `PROJECT_CONTEXT.md` (scNOME-HiC hg38, dropped methylhic_new)

## Latest submitted job (per global CLAUDE.md rule)
`5776408` (scnh_bamproc hg38, pending dependency on 5776407).
Check with `squeue -u jmj7858`.

## Next steps (after 5776408 completes)
1. Verify 188 × `scnomehic/04.bhmem_bam/*.calmd.trinuc_methy.chr21.txt`.
2. Run the scnomehic extract_trinuc line → `summary/trinuc/scnomehic.chr21.txt`.
3. Update `summary/compare_qc_all.py` `load_trinuc_data()` to include
   `methylhic.chr19.txt` and `scnomehic.chr21.txt`.
4. Re-run `compare_qc_all.py` and inspect updated `all_methods_qc.methylation.pdf`.
5. Decide whether to delete `scnomehic/04.bhmem_bam_mm10_wrong/` to reclaim disk.
