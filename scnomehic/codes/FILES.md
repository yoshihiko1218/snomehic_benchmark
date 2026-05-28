# scnomehic/codes — file descriptions

## YAP pipeline (pre-existing)
- `00.rename_fastq_header.sh` — renames raw fastq headers to YAP-expected format
- `01.excute_snakemake.sh` — runs YAP snakemake alignment (bowtie2)
- `01.excute_snakemake_bowtie1.sh` — runs YAP snakemake with bowtie1
- `mapping_config.ini`, `mapping_config_bowtie1.ini` — YAP config files
- `rename_fastq_single.sh`, `rename_fastq_headers_to_yap_format.sh` — fastq renaming helpers

## bhmem pipeline (added 2026-04-18)
- `01.trim.sh` — SLURM array (1-188). trim_galore paired, clip_R1/2=15, three_prime_clip=5,
  -j 32; writes to `03.trimmed_fastq/` and FastQC to `02.fastqc_out/`
- `02.alignment.sh` — SLURM array (1-188). Bhmem alignment on trimmed FASTQ with
  `-nonDirectional -pbat -buffer 100000 -enzymeList dpnII.span_region.bedgraph`;
  writes `.bhmem.bam` + `.flagstat.txt` to `04.bhmem_bam/`
- `03.bamprocess.sh` — SLURM array (1-188). samtools fixmate/markdup/calmd,
  mh_reads_summary.v2.py, Bis-QC (WCH pattern, emits trinuc_methy.chr19/chrM.txt),
  Bis-SNP BisulfiteGenotyper

## Updated Snakemake pipeline rerun (added 2026-05-28)
These replace the standalone bhmem scripts above. They drive the UPDATED
pipeline at `/gpfs/projects/b1198/.../software/sc_NOMeHiC_pipeline` against
workdir `gm_sc_new` (hg38, 188 cells), from scratch through methylation
(no hicluster).
- `dryrun_pipeline.sh` — `snakemake -np` dry-run of the pipeline (resume-aware);
  used to scope what would run before launching.
- `run_pipeline_gm_sc_new.sh` — DRIVER sbatch job (genomicslong, 7-day, 2-core).
  Runs `snakemake --profile slurm --forceall -j 1000`, submitting each rule as a
  child SLURM job (account b1042, partition genomics). Driver log:
  `logs/pipeline_rerun/driver.<jobid>.{out,err}`.

## Logs folder
`logs/{01.trim,02.bhmem,03.bamprocess}/` — SLURM per-task stdout/stderr (old bhmem scripts)
`logs/pipeline_rerun/` — driver job stdout/stderr for the Snakemake rerun
