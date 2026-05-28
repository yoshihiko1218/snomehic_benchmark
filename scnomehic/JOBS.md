# scnomehic SLURM Jobs

## ACTIVE — 2026-05-28: from-scratch rerun via UPDATED Snakemake pipeline
- **Driver job:** `9439299` (`scnh_driver`), partition `genomicslong`, 7-day, 2-core.
  - Submit: `sbatch codes/run_pipeline_gm_sc_new.sh`
  - Driver log: `logs/pipeline_rerun/driver.9439299.{out,err}`
  - Runs `snakemake --profile slurm --forceall -j 1000` on the updated pipeline
    `/gpfs/projects/b1198/.../software/sc_NOMeHiC_pipeline`.
  - **Workdir:** `/gpfs/projects/b1198/.../sc_nomehic_cellline/gm_sc_new` (hg38, 188 cells).
  - **Scope:** trim → align (bhmem) → bamprocess → qc → bisqc → bistools → methylation.
    NO hicluster. `--forceall` regenerates everything (user fixed alignment + bistools).
  - Child jobs: submitted by the driver to account b1042 / partition genomics
    (UUID-named in squeue, ≤48h each). Check with `squeue -u jmj7858`.
  - Pre-flight fix: gunzipped `reference/hg38/...DpnII.span_region.bedgraph(.gz)`
    (mapping rule's `-enzymeList` needs the uncompressed file).
- Status as of submission: PENDING (Priority) on genomicslong.

---

## (historical) earlier standalone-bhmem attempts

## Pipeline outcome
**All in-benchmark bhmem runs produced 0 mapped reads and are INVALID.**
The correct bhmem outputs were already produced externally by the user in:
`/projects/b1198/epifluidlab/yoshii/sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake/`
— 188 cells, each with `.calmd.bam` (~19M mapped reads), `.calmd.trinuc_methy.chr21.txt`,
`.calmd.trinuc_methy.chrM.txt`, etc.

For benchmark QC, we just extract trinuc from that folder — no alignment needed in benchmark/.

## Historical (failed) attempts

### 2026-04-18 — mm10 (wrong genome)
- 5703168 scnh_trim (1-188) COMPLETED (trimmed FASTQs valid, genome-agnostic).
- 5703169 scnh_bhmem mm10 (1-188) COMPLETED, 0 mapped reads (wrong genome).
- 5703170 scnh_bamproc mm10 (1-188) COMPLETED, Bis-QC stopped (empty BAMs).
- Output moved aside to `04.bhmem_bam_mm10_wrong/`.

### 2026-04-19 → 2026-04-22 — hg38 (correct genome, but fastq read names broken)
- 5776407 scnh_bhmem hg38 (1-188) COMPLETED, 0 mapped reads.
- 5776408 scnh_bamproc hg38 (1-188) COMPLETED, Bis-QC stopped (empty BAMs).
- Root cause: `scnomehic/fastq/*-R1.fq.gz` and `*-R2.fq.gz` have trailing `_1`/`_2`
  on read IDs (e.g. `@LH00305:...:1144_1` vs `@LH00305:...:1144_2`). BWA-MEM pairs
  by name up to first whitespace → no pair match → 0 properly paired reads.
- Output in `04.bhmem_bam/` (kept for post-mortem; can delete).

## Root cause
`scnomehic/fastq/*-R{1,2}.fq.gz` were demultiplexed copies with munged read names
(no whitespace-separated field, `_1`/`_2` suffix breaks BWA pairing).
The original raw FASTQs with proper Illumina headers are at
`/projects/b1198/epifluidlab/yoshii/sc_nomehic_cellline/gm_sc_new/00.raw_data/`
and the full processed outputs (demul, trim, bhmem, bamprocess) are in
`gm_sc_new/01.demul_fastq_snakemake/`, `.../03.trimmed_fastq_snakemake/`, and
`.../04.alignment_snakemake/`. Reference pipeline: `/projects/b1198/epifluidlab/yoshii/meningioma_scnomehic/data_batch_1/00-07*.sh`.

## Cleanup candidates (disk reclaim)
- `scnomehic/03.trimmed_fastq/` (~150 GB) — our local trim output, unused
- `scnomehic/04.bhmem_bam/` (hg38 broken BAMs)
- `scnomehic/04.bhmem_bam_mm10_wrong/` (mm10 broken BAMs)
- `scnomehic/02.fastqc_out/` (fastqc from local trim)
