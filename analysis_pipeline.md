# Analysis pipeline — how each QC metric is derived per technology

This documents, **per technology (one folder each)**, how the benchmark's QC metrics
are computed, and the **consistency conventions** that make them biologically
comparable across methods. Metrics: (1) bisulfite conversion, (2) uniquely-mapped
count, (3) MapQ30 mapping rate, (4) per-cell cis-contacts, (5) trans/(trans+cis)
ratio, (6) per-cell cis>1kb contacts, plus HCG / GCH detected-loci.

Master QC notebook: `summary/qc.ipynb` (regenerates all figures in `figures/`).

---

## Consistency conventions (biological definitions held constant across methods)

- **Bisulfite conversion (metric 1):** non-CpG methylation % in the **ACT**
  trinucleotide context (H=A/C/T, non-CpG & non-GpC), reported for **chrM** and an
  **autosome** (chr21 for human / chr19 for mouse). Process may differ (read-level
  trinuc vs allc position-level) — the biological quantity is the same.
- **Fragment = one molecule = one read-name.** Counts collapse R1+R2 to one fragment
  (no double-count). A fragment is uniquely-mapped if ≥1 primary, mapped, non-secondary,
  non-supplementary read carries that name.
- **Duplicates:** removed by **position-based PCR dedup** everywhere
  (`samtools markdup` / `picard MarkDuplicates` / `pairtools dedup` — biologically
  equivalent). Metrics 2 & 3 are reported **both before-dedup and after-dedup**.
- **MapQ30 rate (metric 3) = MapQ30 fragments / mapped fragments** (fraction of mapped
  reads that are high-quality). MapQ threshold = 30.
- **Hi-C contacts (4,5,6):** MapQ≥30, deduplicated, **cis-long threshold = 1 kb**.
  `cis_n = all cis`, `cis_gt1kb = cis with |pos1-pos2|>1000`,
  `trans_ratio = trans / (cis + trans)`.
- **HCG/GCH loci:** count of **≥1×-covered** cytosine sites, **destranded** (collapse
  ± to the + strand CpG/GpC), **GCG removed** by reference lookup. Caller may differ.

### Consistent recompute of metrics 2 & 3 (all technologies)
`summary/frag_counts.py` (per-cell) + `summary/frag_jobs/frag_array.sh` (SLURM array,
markdup-if-needed) → `summary/frag_counts/<ds>/<cell>.tsv`, collected by
`summary/frag_jobs/collect.py` → **`summary/frag_counts_all.tsv`** with columns
`uniq_preDedup, mapq30_preDedup, rate_preDedup, uniq_postDedup, mapq30_postDedup, rate_postDedup`.
BAM used per technology is listed in each section below and in `summary/frag_jobs/JOBS.md`.

---

## nagano/  — scHi-C (Nagano 2013, mm10, 15 cells, pure Hi-C)
Pipeline: fastp → **bwa-mem2 `-SP5M`** → `samtools markdup` → `mh_reads_summary.v2.py`.
- **Metrics 2,3:** `frag_counts.py` on `alignment/{cell}.markdup.bam` (dup-flagged).
- **Metrics 4,5,6:** parse `alignment/{cell}.summary.txt` fields
  `UniqMappedMapQ30NoPcrCis` (cis), `UniqMappedMapQ30NoPcrCisMore1kb` (cis>1kb),
  `UniqMappedMapQ30NoPcrTrans` (trans). trans_ratio = trans/(cis+trans).
- No methylation (metrics 1 / loci N/A). Note: high trans_ratio (~0.17) is a real
  quality signature of this early scHi-C assay.

## droplethic/  — Droplet Hi-C (Chang 2025, hg38, ~688k barcodes; 3,668 valid)
Pipeline: Rupture (Trim Galore → bowtie1 barcode → **bwa mem `-SP5M`** → CB-tag →
**pairtools** parse/sort/dedup). Per-barcode QC by `my_project/rupture_qc.py`.
- **Metrics 2,3:** `frag_droplethic.py` streams the CB-tagged BAM
  `my_project/03.mapping/SRR27586278_hg38.bam`, per-barcode.
- **Metrics 4,5,6:** `my_project/SRR27586278_hg38.per_cell_qc.valid.tsv`
  (`UniqMappedMapQ30NoPcrCis`/`...Cis1kb`/`...Trans`).
- ⚠️ **DEDUP CAVEAT:** the `pairtools dedup` step was not run (`pairparse.txt`
  `total_dups=0`); the CB-tagged BAM has no duplicate flags, so droplethic values are
  currently **pre-dedup** (before==after). To make it fully comparable, run
  `summary/frag_jobs/droplethic_dedup.sbatch` (per-barcode dedup over the 1.17 TB
  position-sorted pairsam; same position-based definition as the others). Kept as-is
  by decision 2026-06-30.
- No methylation.

## scnomehic/  — scNOMe-HiC (this method; GM12878, hg38, 188 cells → 187 QC-pass)
Pipeline: Trim Galore → **bhmem** (bisulfitehic BWA-MEM) → `samtools markdup`/calmd →
`mh_reads_summary.v2.py` + Bis-QC/Bis-SNP. **Real bhmem output is EXTERNAL** at
`/home/jmj7858/epifluidlab/sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake/`
(local `scnomehic/alignment/` is the YAP bowtie2 comparison arm).
- **Metrics 2,3:** `frag_counts.py` on external `gm_sc_new/*.calmd.bam` (dup-flagged);
  QC-passing cell set = `.../s1/QC/gm_passed.txt` (187).
- **Metrics 4,5,6:** external `gm_sc_new/{cell}.summary.txt.gz`
  (`UniqMappedMapQ30NoPcrCis`/`...CisMore1kb`/`...Trans`).
- **Metric 1:** external `s1/QC/chrM_gch_hcg/gm.{chrM,chr21}.txt.gz` (`noncpg`).
- **Loci:** BisSNP NOMe HCG/GCH `6plus2` → `summary/gch_hcg_counts/scnomehic.{hcg,gch}`.

## snmCseq3/  — snm3C-seq (Liu 2023, mm10, 100 cells → 98 QC)
Pipeline: cutadapt → **two aligners** (YAP bismark+bowtie2 m3C `alignment/`, and
bhmem `04.bhmem_bam/`). Both consumed by qc.ipynb.
- **Metrics 2,3:** `frag_counts.py` on `04.bhmem_bam/{cell}.calmd.bam` (bhmem,
  markdup-flagged). (2 cells excluded: empty calmd.bam.)
- **Metrics 4,5,6 — TWO MapQ versions kept:**
  - **MapQ≥10 (YAP native):** `alignment/stats/MappingSummary.csv.gz`
    `CisShortContact`+`CisLongContact` (cis), `CisLongContact` (cis>1kb, YAP CisLong
    = 1 kb via `min_gap=1000`), `TransContact`.
  - **MapQ≥30 (matched to other Hi-C):** re-filter `alignment/Group*/bam/{cell}.3C.sorted.bam`
    at `samtools view -q 30`, re-run `yap-internal generate-contacts --min_gap 1000`
    (`mapping` conda env) → `.counts.txt` (`CisShort,CisLong,Trans`). Job in
    `summary/frag_jobs/snmCseq3_contacts_q30.sh` →
    `summary/frag_counts/snmCseq3_contacts_q30/{cell}.q30.counts.txt`.
  - cis_n=CisShort+CisLong; cis_gt1kb=CisLong; trans_ratio=Trans/(CisShort+CisLong+Trans).
  - NB: m3C contacts are chimeric split-read (multi-way) — a method-inherent difference
    from the paired-mate contacts of nagano/scnomehic/droplethic.
- **Metric 1:** `summary/trinuc/snmCseq3.{chrM,chr21}.txt` (from bhmem
  `04.bhmem_bam/*.calmd.trinuc_methy.*` via `summary/extract_trinuc.py`).
- **Loci:** allcools allc (`alignment/stats/AllcPaths.tsv`) → destrand + GCG removal →
  `summary/gch_hcg_counts/snmCseq3.hcg`.

## scnome/  — scNOMe-seq (Pott 2017, hg38, 23 cells: 12 GM12878 + 11 merged K562)
Pipeline: Trim Galore → **Bismark SE per-mate** (`--non_directional`) → markdup →
`coverage2cytosine --nome-seq`.
- **Metrics 2,3:** `frag_counts.py` on `04.alignment/{cell}_{1,2}.rmdup.RG.bam`
  (per-mate R1+R2, dup-flagged; collapse to fragments by read-name).
- **Metric 1:** `scnome_qc_summary.csv` `chrM_noncpg`, `chr21_noncpg` (per-cell
  Combined; Bismark trinuc).
- **Loci (HCG + GCH; NOMe):** per-cell **union** of the two mates' destranded loci,
  `summary/scnome_loci_percell.py` → `summary/gch_hcg_counts/scnome.{hcg,gch}`.
- No Hi-C.

## snmCseq2/  — snmC-seq2 (Luo 2018, hg38+mm10, 249 cells; **mm10 subset = 96** used)
Pipeline: cutadapt → **Bismark SE per-mate** → dedup → methyl-extractor. Also a **YAP**
`mc` re-run (`yap_mapping_{hg38,mm10}`).
- **Metrics 2,3:** `frag_counts.py` on raw `05.align/{cell}_{1,2}.clean_bismark_bt2.bam`
  + **markdup** (so both dedup versions are available; the YAP `final.bam` has dups
  already removed → after-dedup only, which is why raw+markdup is used here). mm10 cells
  (`codes/cells_mm10.txt`).
- **Metric 1:** `chrM` from `snmcseq2_qc_summary.csv` `chrM_noncpg`; **autosome (chr19)**
  from `summary/trinuc/snmCseq2.chr19.txt` (per-mate, **averaged per cell**).
- **Loci (HCG):** Bismark cov `06.methy/{cell}.clean_bismark_bt2.deduplicated.bismark.cov.gz`
  (mm10) → destrand + GCG removal.
- No Hi-C / no GCH (not NOMe).

## smallwood/  — scWGBS (Smallwood 2014, mm10, 51 cells)
Pipeline: Trim Galore → **Bismark** (hg38 human-contamination depletion → mm10 SE) →
markdup → cov + BisSNP.
- **Metrics 2,3:** `frag_counts.py` on `05.align_mm10/{cell}.rmdup.RG.bam` (dup-flagged, SE).
- **Metric 1:** `smallwood_qc_summary.csv` `chrM_noncpg`, `chr19_noncpg`.
- **Loci (HCG):** Bismark cov `06.methy/{cell}.dedup.bismark.cov.gz`.
- No Hi-C / no GCH (not NOMe).

## snmCAT/  — snmCAT-seq / snmC2T (Luo 2022, human brain, hg38, 100 cells → 99)
Pipeline: YAP `mct --nome` (**bismark** methylome + **STAR** transcriptome + allcools).
- **Metrics 2,3:** `frag_counts.py` on `mapping_brain/Group*/bam/{cell}.dna_reads.bam`
  + **markdup** (dna_reads.bam has no dup flags).
- **Metric 1:** computed from YAP **allc** via `summary/extract_trinuc_snmCAT.py`
  (ACT/ACG/GCT % per chrom, tabix) → `summary/trinuc/snmCAT.{chrM,chr21}.txt`.
  NB: autosomal non-CpG is elevated by **real brain neuronal mCH**, not conversion failure.
- **Loci (HCG + GCH; NOMe):** allcools allc (`mapping_brain/stats/AllcPaths.tsv`) →
  destrand + GCG removal → `summary/gch_hcg_counts/snmCAT.{hcg,gch}`.
- No Hi-C.

## methylhic/ , methylhic_new/  — excluded from the current benchmark figures (per user).

---

## Aggregators / scripts (in `summary/`)
- `collect_conversion.py` → `conversion_percell.csv` (metric 1, per-cell, never per-mate).
- `extract_trinuc.py`, `extract_trinuc_snmCAT.py` → `summary/trinuc/*` (conversion inputs).
- `frag_counts.py` + `frag_jobs/` (manifests, `frag_array.sh`, `frag_droplethic.py`,
  `collect.py`) → `frag_counts_all.tsv` (metrics 2,3, both dedup).
- `snmCseq3_contacts_q30.sh` → snmCseq3 MapQ30 contacts.
- `make_dataset_summaries.py` → `gch_hcg_counts/all_methods.summary.txt` (loci).
- `dist_check.py` → cross-technology distribution sanity check.
- `qc.ipynb` → all figures.

## Known method-inherent caveats (not fixable by reprocessing)
1. Hi-C contact **model** differs: snm3C-seq = chimeric split-read (multi-way);
   nagano/scnomehic/droplethic = paired-mate. Thresholds (1kb/MapQ30/dedup) are matched.
2. **snmCAT** autosomal non-CpG conversion is confounded by real brain mCH.
3. Conversion autosome proxy is chr21 (human) vs chr19 (mouse).
4. **droplethic** not yet PCR-deduplicated (see its section).
