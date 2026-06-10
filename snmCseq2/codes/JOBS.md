# JOBS — snmCseq2

## 05.methy_extract — Bismark CpG methylation extraction (snmC-seq2)
- **Job name:** sc2methy
- **Job ID:** 4287414  (array 1-498%60, one task per `_1`/`_2` mate-sample in
  `methy_samples.txt`; the _1/_2 split is kept, matching the 498 BisSNP samples)
- **Submitted:** 2026-06-10
- **Command:** `sbatch codes/05.methy_extract.sh`
- **Logs:** `snmCseq2/logs/05.methy_extract/sc2methy.<arrayid>.txt` (+ `.err`)
- **Genome:** hg38 (`hg38_bismark/`); snmC-seq2 is NOT NOMe -> CpG only, no GCH.
- **What it does:** from each `05.align/<sample>.clean_bismark_bt2.bam` (raw
  bismark alignment; the only bam present for all 498):
  - Step 0: `deduplicate_bismark -s` (protocol: dedup before extraction; uniform
    across all 498 since some rmdup bams were cleaned up).
  - Step 1: `bismark_methylation_extractor -s --comprehensive --bedGraph
    --genome_folder hg38_bismark` -> `06.methy/<sample>.clean_bismark_bt2.deduplicated.bismark.cov.gz`
    (one row per covered CpG; per-locus, per-strand, no SNP filter — YAP/Bismark
    convention, NOT BisSNP). Detected CpG (HCG) loci = #rows.
  - Temp `*.deduplicated.bam` removed after the cov is written.
- **Validation:** SRR6911624_1 -> 956,051 CpG loci (deep cell); pipeline ran
  dedup + extractor cleanly.
- **Note:** snmCseq2 was aligned with bismark `-pbat`; XM tags are correct, `-s`
  extraction is appropriate.
- **Next:** count loci into `summary/gch_hcg_counts/` after completion.

### Status checks
- `squeue -u jmj7858 | grep sc2methy`
- On failure inspect `logs/05.methy_extract/sc2methy.<id>.err`.

## 07.yap_mc — yap (cemba_data) mc mapping, MIXED-SPECIES (snmC-seq2)
- **Why:** reprocess snmCseq2 with `yap` like snmCAT (mct) / snmCseq3 (m3c).
  Mode = `mc` (methylation-only). yap 1.6.9, conda env `mapping`.
- **KEY FACT:** snmCseq2 is a MIXED-SPECIES dataset — 153 hg38 cells + 96 mm10 cells
  (per `snmCseq2_genome_map.tsv`; each cell's _1/_2 share one genome). yap takes one
  bismark ref per run, so there are TWO independent yap runs.
- **Job IDs:**
  - `sc2_yap_hg38` = **4327550** (array 1-64, output `yap_mapping_hg38/`)
  - `sc2_yap_mm10` = **4327551** (array 1-64, output `yap_mapping_mm10/`)
- **Submitted:** 2026-06-10
- **Submit cmds:**
  - `sbatch --job-name=sc2_yap_hg38 --export=ALL,GENOME=hg38 codes/07.run_yap_snakemake.sh`
  - `sbatch --job-name=sc2_yap_mm10 --export=ALL,GENOME=mm10 codes/07.run_yap_snakemake.sh`
- **Setup steps (already done):**
  1. `bash codes/06.yap_symlink.sh` -> `fastq_yap_{hg38,mm10}/*-R[12].fq.gz` + `codes/cells_{hg38,mm10}.txt`
  2. `yap default-mapping-config --mode mc -v V2 ...` -> `codes/mapping_config_yap_{hg38,mm10}.ini`
  3. `yap start-from-cell-fastq --output_dir yap_mapping_<g> --config_path ... --fastq_pattern "fastq_yap_<g>/*-R[12].fq.gz"`
     -> 64 Group Snakefiles each + `snakemake/snakemake_cmd.txt`
- **Refs:** hg38_bismark + GCA_000001405.15_GRCh38_no_alt fasta/chrom.sizes;
  mm10_bismark + mm10.fa/mm10.chrom.sizes.
- **Logs:** `codes/logs/07.yap_mc/sc2_yap_<genome>.<arrayid>.{out,err}`
- **After both arrays finish:** `yap summary yap_mapping_hg38/` and `yap summary yap_mapping_mm10/`
  -> `MappingSummary.csv.gz` in each run dir.
- **NOTE:** an earlier single-genome (mm10-only) attempt (job 4327205) was CANCELLED once
  the mixed-species map was discovered; its artifacts are in `.old_single_genome_attempt/`.

### Status checks
- `squeue -u jmj7858 | grep sc2_yap`
- On failure: `codes/logs/07.yap_mc/sc2_yap_<genome>.<id>.err`

### 07.yap_mc — FIX 1 (bismark_reference NameError) + resubmit
- **Bug:** first run (4327550/4327551) finished in ~3 min producing ZERO bam/allc.
  Cause: `cemba_data/mapping/pipelines/mc.py` drops `bismark_reference` from the
  Snakefile header when the config still contains a `hisat3n_dna_reference` key
  (yap `default-mapping-config` leaves a placeholder line). The mc.Snakefile template
  calls `bismark {bismark_reference}` -> `NameError: name 'bismark_reference' is unknown`.
- **Fix:** `codes/08.fix_bismark_config.sh` deletes the `hisat3n_dna_reference` line from
  all 4 configs (2 source + 2 copied), then `yap update-snakemake -o yap_mapping_<g>`
  regenerates Snakefiles. Verified header now has `bismark_reference = '<path>'` (line 20).
- **Resubmitted:** sc2_yap_hg38 = **4330289**, sc2_yap_mm10 = **4330290** (array 1-64).
