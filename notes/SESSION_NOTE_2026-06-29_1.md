# Session Note 2026-06-29 (1)

## Goal
Clean each benchmark "by how it is ran": go dataset-by-dataset, identify the
canonical run path vs dead-end experiments / stale intermediates, document the
real pipeline, and produce **delete commands for the user to run by hand**
(I have no delete permission — only propose `rm` commands after discussion).

User decisions (this session):
- Clean = remove dead-end experiments + document real run path + free disk space.
- Do NOT standardize folder layout (Tier-2, path-coupling risk).
- Disposal = discuss, then hand the user delete commands to run manually.
- Pace = one benchmark at a time, user reviews each.

Proposed order: snmCseq3 → snmCseq2 → snmCAT → scnome → smallwood → methylhic
→ methylhic_new → scnomehic → nagano → droplethic.

## Progress log
- (start) Read README.md, REORG_PLAN.md, PROJECT_CONTEXT.md for context.

### Benchmark 1/10: snmCseq3 — DONE (proposed, awaiting user run)
- Verified via repo-wide grep that qc.ipynb consumes BOTH `alignment/` (YAP,
  MappingSummary) AND `04.bhmem_bam/` (bhmem summaries/trinuc/CG beds). Both KEEP.
  (Trap avoided: PROJECT_CONTEXT said "YAP only" but bhmem is also consumed.)
- Dead-end experiment dirs (bhmem-vs-YAP NM/XG/MAPQ investigation) are
  unreferenced by any QC pipeline → safe to delete (~275 GB).
- Wrote `snmCseq3/FILES.md` documenting canonical run vs dead-ends.
- Wrote `cleanup/cleanup_snmCseq3.sh` (delete script for user to run by hand).
- TODO after user deletes: drop 5 `alignment_bowtie1/_mapq0` git-add lines from
  `tools/commit_all.sh` (lines ~24-33).

### Pipeline survey (all 10) — DONE via 10 parallel Explore agents
Established the actual tool chain + canonical output + dead-end arms per dataset.
Key cross-cutting facts:
- Three pipeline families: (A) YAP/cemba_data (bismark+bowtie2+allcools; m3c/mc/mct)
  → alignment/stats/MappingSummary.csv.gz; (B) Bismark-SE per-mate methylation
  (trim_galore/cutadapt + bismark SE + methyl-extractor/BisSNP); (C) pure Hi-C
  (bwa-mem2/bwa -SP5M + pairtools + mh_reads_summary.v2.py).
- bhmem = bisulfitehic BWA-MEM aligner = "the method"; used as method aligner
  (scnomehic, real output EXTERNAL at b1198/.../gm_sc_new/04.alignment_snakemake)
  and as a comparison arm in methylhic (04.alignment), snmCseq3 (04.bhmem_bam),
  methylhic_new (legacy codes).
- MANY datasets ran 2 pipelines (canonical + comparison arm). Must not assume the
  numbered/standard dir is canonical — verify against qc.ipynb each time.
- scnomehic: local alignment/ (YAP bowtie2) canonical; alignment_bowtie1/ (341G)
  is a DEAD-END param experiment; local 04.bhmem_bam = INVALID (0 mapped, broken
  read names) — real bhmem is external.
- snmCAT: original H1/HEK snmCT download already removed; current fastq_brain/
  mapping_brain IS genuine brain snmC2T-seq NOMe, IS used in qc.ipynb.

### Benchmark 2/10: scnomehic — DONE (proposed, awaiting user run)
- Verified: qc.ipynb consumes `alignment/stats/MappingSummary.csv.gz` (local YAP
  bowtie2). `alignment_bowtie1/` (341G) referenced 0x → dead-end.
- bhmem (method) valid output is EXTERNAL on b1198 (gm_sc_new) — untouched.
- Wrote `scnomehic/FILES.md` + `cleanup/cleanup_scnomehic.sh` (frees ~341 GB).

### Summary-page (qc.ipynb) provenance audit — DONE
Traced every dataset → pipeline → file consumed by qc.ipynb. Key results:
- CORRECTION to snmCseq3: cells 32/33 OVERWRITE snmCseq3 with YAP MappingSummary
  before plots; so final figures use YAP for align+contacts (metrics 2-6), bhmem
  only for conversion (1) + HCG-loci beds. Fixed snmCseq3/FILES.md table.
- Figure→dataset sets: base-align=7 (nagano,droplethic,scnome,smallwood,snmCseq2,
  snmCseq3,scnomehic); Hi-C contacts=4 (nagano,droplethic,snmCseq3,scnomehic);
  loci HCG=5 methylation, GCH=NOMe only (scnome,snmCAT,scnomehic).
- methylhic + methylhic_new are NOT in current qc.ipynb (only older compare_qc_*).
  Flag before cleaning those two.
- scnomehic figure data is ENTIRELY external (b1198 gm_sc_new); local folder only
  contributes via the (kept) alignment/ — reinforces alignment_bowtie1 = dead weight.
- cell 28 snmCseq3 conversion; cell 5 nagano align summary; cell 8 droplethic
  valid.tsv; cell 49 hic_df (4 datasets); cell 57 gch_hcg loci from
  summary/gch_hcg_counts/all_methods.summary.txt (built by make_dataset_summaries.py).

### User directives (2026-06-30)
- IGNORE methylhic + methylhic_new for cleanup (out of scope).
- User confirms they ran snmCseq2 with YAP too → yap_mapping_hg38/mm10 are
  deliberate; KEEP (not dead-ends).

### Benchmark 3/10: snmCseq2 — DONE (proposed, awaiting user run)
- Two arms: Bismark-SE (canonical → snmcseq2_qc_summary.csv) + YAP (yap_mapping_*,
  kept). Provenance: conversion=trinuc(05.align); count/%=qc_summary(Bismark);
  HCG loci=06.methy cov.gz; HCG sanity=yap_mapping_mm10 AllcPaths.
- HEADLINE: 06.methy/{CpG,CHG,CHH}_context_*.txt = 612G of bismark intermediates,
  unreferenced, cov.gz/trinuc/qc_summary already exist → delete frees ~612G.
- Wrote snmCseq2/FILES.md + cleanup/cleanup_snmCseq2.sh.
- 249 cells (153 hg38 + 96 mm10).

### Conversion-figure collection + plotting (2026-06-30) — data layer DONE
Goal: add the missing metric-1 (bisulfite conversion) cross-dataset figure, 6
datasets incl snmCAT, chrM + chr21/chr19; also add snmCAT to metrics 2&3.
New reproducible scripts in summary/:
- extract_trinuc_snmCAT.py -> trinuc/snmCAT.{chrM,chr21}.txt (ACT/ACG/GCT % from
  YAP allc via tabix; 4-mer context, trinuc = first 3 chars; pos-level Σmc/Σcov).
- snmCseq3.chrM.txt via existing extract_trinuc.py on 04.bhmem_bam (98 cells).
- collect_conversion.py -> conversion_percell.csv (554 cells). Per-CELL never
  per-mate. Cell sets: scnome23, smallwood51, snmCseq2 96(mm10), snmCseq3 98,
  scnomehic 187(passed), snmCAT 99. snmCseq2 chr19 from per-mate trinuc averaged;
  chrM from qc_summary. scnomehic from external gm.{chrM,chr21}.
- plot_conversion.py -> figures/conversion_{chrM,chr21chr19}_violin.{pdf,png}.
Medians non-CpG%: chrM scnome6.25/scnomehic2.05/snmCAT1.58/smallwood0.36/
  snmCseq2~0/snmCseq3~0; auto snmCAT3.35/snmCseq3 1.69/snmCseq2 1.63/scnomehic1.10/
  smallwood0.80/scnome0.35. (snmCAT high = real brain mCH.)
TODO (needs qc.ipynb edit + rerun): add snmCAT to readcount(cell48)+mapq30(cell47)
  using mean(R1,R2 UniqueMappedReads) & mean(R1,R2 MappingRate) [no MapQ30 cols in
  snmCAT MappingSummary]; optionally fold conversion cells into qc.ipynb.

### Conversion-figure + snmCAT integration — DONE (2026-06-30)
- Edited qc.ipynb (programmatic JSON edit; backup at
  summary/_backups/qc.ipynb.pre_snmCAT_conversion.bak): inserted snmCAT-build cell
  after df_all_gm; added snmCAT to dataset_order/palette/concat in mapq30 + readcount
  cells; appended 2 conversion-figure cells.
- Executed qc.ipynb end-to-end (jupyter nbconvert --execute --inplace) -> rc=0.
  Regenerated all figures incl new figures/conversion_{chrM,chr21chr19}_violin.*
  and snmCAT now in mapq30 (~75%) + readcount panels (8 datasets).
- snmCAT mapped N=99 cells. Conversion n per dataset verified.
- Documented in summary/CONVERSION_README.md.
- NEXT: resume cleanup at benchmark 4/10 = scnome (per earlier plan), and/or commit.
