# FILES — summary/

QC summary, cross-technology comparison, and the `qc.ipynb` notebook.

## Benchmark data summary workbook
- **build_data_summary_xlsx.py** — builds `benchmark_data_summary.xlsx` from the
  facts distilled across every dataset folder's `FILES.md`/`RESULTS.md`, the
  top-level `README.md`, `PROJECT_CONTEXT.md`, and `analysis_pipeline.md`.
- **benchmark_data_summary.xlsx** — Excel summary of the data used in the benchmark.
  Sheets: **Datasets** (one row per dataset: method, publication, **publication DOI**,
  genome, cell type, modalities, cell counts, **GEO accession**, SRA accessions,
  pipeline, notes), **QC metrics** (the 6 metrics + HCG/GCH loci), **Caveats**
  (method-inherent caveats). GEO/DOI verified against NCBI GEO/BioProject and the
  source publications (2026-07).

## Notebook
- **qc.ipynb** — master QC notebook. Loads per-cell metrics for every method
  (nagano, droplethic, scnome, smallwood, snmCseq2, snmCseq3, scnomehic) and
  draws the publication violin panels. Recent additions:
  - cis > 1 kb rate violin (Hi-C block).
  - **GCH / HCG detected-loci** violins (see below).

## GCH / HCG detected-loci counting (added this work)
"Detected loci" = number of data rows in a BisSNP `*.6plus2.bed` file (one row
per detected cytosine in that context; minus the single `track` header line).
All datasets use the **same caller (BisSNP-0.90)** and all 6plus2 files are
**destranded** (`+` strand only), so counts are directly comparable.

- **count_gch_hcg.py** — counts GCH and HCG loci per cell across all methylation
  methods, writing to `gch_hcg_counts/`. Per-dataset sources:
  - scnomehic: `gm_sc_new/07.bistools_snakemake/methylation/<cell>/<cell>.cyt.filtered.sort.{GCH,HCG}.6plus2.bed.gz` (GM12878, **BisSNP-1.0.1** — the newer caller; per-strand, ~2.25x the older BisSNP-0.90 `calmd.cytosine.filtered.sort` output)
  - scnome:    `scnome/04.alignment/<cell>.rmdup.RG.cytosine.filtered.sort.{GCH,HCG}.6plus2.bed`
  - smallwood: `smallwood/05.align_mm10/<cell>.rmdup.RG.cpg.filtered.sort.CG.6plus2.bed` (CG → HCG; no GCH)
  - snmCseq2:  `snmCseq2/05.align/<cell>.rmdup.RG.cpg.filtered.sort.CG.6plus2.bed` (CG → HCG; _1/_2 split)
  - snmCseq3:  `snmCseq3/04.bhmem_bam/<cell>.calmd.cpg.filtered.sort.CG.6plus2.bed` (CG → HCG; no GCH)
  GCH (GpC accessibility) exists only for the NOMe methods (scnome, scnomehic).
- **count_gch_hcg.log** — run log.
- **gch_hcg_counts/all_methods.gch_hcg_loci.txt** — combined per-cell table:
  `dataset, sample, GCH_n, HCG_n`.
- **gch_hcg_counts/<dataset>.gchhcg.txt** — per-dataset `sample, GCH_n, HCG_n`.

### snmCseq3 cross-check vs YAP/allcools (BisSNP vs YAP)
- **count_yap_cpg.sh** — counts CpG loci per snmCseq3 cell from the YAP ALLC
  files (`snmCseq3/alignment/.../allc/<cell>.allc.tsv.gz`, context starts `CG`).
- **count_yap_cpg.log** — run log.
- **gch_hcg_counts/snmCseq3_yap_cpg.tsv** — YAP CpG counts per cell.
- **gch_hcg_counts/snmCseq3_yap_vs_bissnp.tsv** — merged comparison (65 cells).
  Result: Pearson **r = 0.9995**; YAP is **4.26×** BisSNP, decomposed as
  ~2× (YAP per-strand vs BisSNP destranded) × **2.13×** (BisSNP SNP/quality
  filtering keeps ~47% of covered CpGs). Conclusion: caller choice does not
  change relative comparisons; keep all datasets on BisSNP for the figure.

## HCG per-technology calculation (tracked in each tech's codes/)
HCG = covered CpG with **GCG removed** (= ACG/CCG/TCG), one definition for all.
How GCG is removed differs by how BisSNP was run, and each technology has its own
script so the calculation is tracked:

- **summary/hcg_lib.py** — shared core: genome load, `classify()` (HCG/GCG from
  reference), `count_6plus2()` (manual GCG removal), `count_rows()` (NOMe files
  that already exclude GCG).
- **summary/count_hcg_destranded.py** — the consistent engine: per-cell, collapse
  +/- strand observations of each CpG (destrand) and remove GCG by reference-genome
  context. `compute(dataset)` is called by each tech's compute_hcg.py.
- **<tech>/codes/compute_hcg.py** — per technology; documents aligner + caller and
  writes `summary/gch_hcg_counts/<tech>.hcg.txt`. CONSISTENT calling + GCG removal
  for the 4 non-scnomehic datasets:
  - scnome    : **Bismark** methylation extraction -> destrand + genome GCG removal (hg38)
  - smallwood : **Bismark** methylation extraction -> destrand + genome GCG removal (mm10)
  - snmCseq2  : **Bismark** (-pbat) extraction -> destrand + genome GCG removal (hg38/mm10)
  - snmCseq3  : **allcools** (bhmem; no Bismark) -> destrand + genome GCG removal (mm10)
  - scnomehic : **BisSNP-1.0.1 NOMe** HCG.6plus2 as-is (the exception; summary/compute_hcg_scnomehic.py)
- **summary/submit_compute_hcg.sh** — SLURM array (1-4) for the consistent ones.
- **summary/assemble_hcg_all.py** — combines the 5 `<tech>.hcg.txt` ->
  `gch_hcg_counts/all_methods.hcg_loci.txt` (dataset, sample, HCG_n).
- Median HCG: scnome 860,392 | smallwood 2,637,415 | snmCseq2 735,594 |
  snmCseq3 1,165,000 | scnomehic 1,241,579.

WHY Bismark/allcools (not BisSNP) for the 4: BisSNP applies a -mmq 30 mapping-
quality cut that drops ~95% of low-MAPQ PBAT reads (snmCseq2 _2 mates) very
unevenly (~3-21x fewer loci than Bismark; see comparison tsvs). Bismark/allcools
have no such cut, so calling AND GCG removal are consistent across the 4.
NOTE: scnomehic (BisSNP NOMe) still carries the MAPQ>=30 cut -- it's the one
exception, kept on its native NOMe pipeline per user.

## Pre-existing comparison artifacts (not generated here)
- compare_qc*.py / *.csv / *.pdf, all_methods_qc.*, extract_trinuc.py,
  run_extract_trinuc.sh, trinuc/ — earlier QC pipeline outputs.

## Subfolders added during reorg (2026-06-18)
- **_archive/** — one-off notebook-patcher scripts (`_add_*.py`, `_fix_*.py`,
  `_update_*.py`, `_ylim0_nb.py`). Superseded; kept for reference, safe to delete.
- **_backups/** — `qc.ipynb.bak_*` snapshots taken before edits (gitignored).
- **trinuc/**, **hic_cache/**, **gch_hcg_counts/** — real inputs read by `qc.ipynb`
  (do not move; the notebook references these relative paths).
