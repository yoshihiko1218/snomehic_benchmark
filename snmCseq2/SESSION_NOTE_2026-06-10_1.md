# Session Note 2026-06-10 (1) — Run snmC-seq2 through yap (mc mode)

## Goal
Re-process snmCseq2 (snmC-seq2, mouse mm10, methylation-only, 249 cells) with `yap`
(cemba_data), mirroring how snmCAT (mct) and snmCseq3 (m3c) were run. Mode = `mc`.

## Context discovered
- yap 1.6.9 in conda env `mapping` (/projects/b1198/.../conda/envs/mapping/bin/yap)
- Existing snmCseq2 pipeline used bismark directly (codes/02.alignment.sh -> 05.align, 06.methy).
  yap run is NON-DESTRUCTIVE: new folder `yap_mapping/`, fastq symlinks in `fastq_yap/`.
- 249 paired cells: 01.fastq/SRRxxxxxxx_{1,2}.fq.gz
- mm10 refs: /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark (Bisulfite_Genome + mm10.fa),
  mm10/mm10.fa, mm10/mm10.chrom.sizes
- Proven recipe (from snmCAT/codes/01,02 and snmCseq3/codes/01):
  symlink *-R[12].fq.gz -> default-mapping-config -> start-from-cell-fastq -> SLURM array of
  snakemake_cmd.txt lines (with `--rerun-incomplete --unlock` then `--rerun-incomplete`) -> yap summary.

## Actions (this session)
- Explored snmCAT (mct) + snmCseq3 (m3c) yap recipes; identified mode `mc` for snmC-seq2.
- Wrote codes/06.yap_symlink.sh, codes/07.run_yap_snakemake.sh.
- **DISCOVERY (critical):** snmCseq2 is MIXED-SPECIES. snmCseq2_genome_map.tsv =
  153 hg38 cells + 96 mm10 cells. Existing bismark bam header (chr1 LN:248956422)
  confirmed hg38 for SRR6911624. The working-copy 02.alignment.sh (mm10) was misleading.
- First submitted a single mm10-only array (job 4327205) -> CANCELLED on discovering
  the mix. Artifacts moved to .old_single_genome_attempt/.
- Restructured into TWO yap runs:
  - bash codes/06.yap_symlink.sh -> fastq_yap_hg38 (153 cells/306 links),
    fastq_yap_mm10 (96 cells/192 links), codes/cells_{hg38,mm10}.txt
  - yap default-mapping-config (mc, V2) -> codes/mapping_config_yap_{hg38,mm10}.ini
  - yap start-from-cell-fastq -> yap_mapping_{hg38,mm10}/ (64 Group Snakefiles each)
- SUBMITTED: sc2_yap_hg38 = **4327550**, sc2_yap_mm10 = **4327551** (both array 1-64).
- Recorded in codes/JOBS.md. Now monitoring (frequent at first, backing off when stable).


## 18:54 update — queued, waiting on cluster
- FIX applied & verified: bismark_reference NameError (see JOBS.md FIX 1). Snakefile header
  now defines `bismark_reference`. Resubmitted hg38=4330289, mm10=4330290.
- Genomics partition 100% full (149A/0I/2O). Est start ~22:27 (hg38) / 22:34 (mm10).
- Nothing to debug; pure queue wait. Monitoring cadence widened to ~1h, converging near 22:27.
