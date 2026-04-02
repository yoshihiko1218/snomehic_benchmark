# scNOME-HiC Benchmark — Project Context

> **Auto-maintained by AI agent. Last updated: 2026-03-26**
> New chat: paste the retrieval command at the bottom of this file.

---

## Project Goal

Benchmark comparison of single-cell HiC and/or methylome sequencing methods.
Produce QC metrics that share a common base (alignment) + modality-specific matrices (methylation, HiC contacts).

**Cluster path:** `/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/`
**Local path:** `/home/jmj7858/epifluidlab/workspace/scnomehic_paper/benchmark/`
**Conda env:** `scnomehic`
**SLURM account:** `b1042`, partition `genomics`

---

## Methods in Benchmark (9 total)

| Folder | Label | Pipeline type | Genome | Notes |
|---|---|---|---|---|
| `methylhic/` | MethylHiC | YAP (m3C) | mm10 | 59 cells |
| `methylhic_new/` | MethylHiC-new | YAP (m3C) | mm10 | 96 cells |
| `snmCseq3/` | snmC-seq3 | YAP | mm10 | 100 cells |
| `scnomehic/` | scNOME-HiC | YAP | mm10 | 376 cells |
| `smallwood/` | Smallwood | Trim Galore + Bismark SE | mm10 | 51 cells |
| `scnome/` | scNOMe | Trim Galore + Bismark SE (R1+R2) | hg38 | 41 cells |
| `snmCseq2/` | snmC-seq2 | cutadapt + Bismark SE (R1+R2) | hg38+mm10 | 249 cells (153 hg38, 96 mm10) |
| `nagano/` | Nagano | fastp + Bowtie2 PE | mm10 | 15 cells |
| `droplethic/` | DropletHiC | CB-tagged BAM + pairs | hg38 | ~688k barcodes |

---

## QC Infrastructure

### Per-method QC scripts

**YAP methods** (methylhic, methylhic_new, snmCseq3, scnomehic):
- Input: `alignment/stats/MappingSummary.csv.gz` — generated automatically by YAP

**Smallwood** (`smallwood/`):
- `smallwood_qc_per_cell.py` — parses Trim Galore, Bismark SE, BAM stats, trinuc, CpG bed, BisQC
- `collect_smallwood_qc.py` — aggregates to `smallwood_qc_summary.csv`
- `run_qc.sh` — SLURM array job (one task per cell), then collect

**scNOMe** (`scnome/codes/`):
- `scnome_qc_per_cell.py` — parses Trim Galore, Bismark SE, BAM summary, 6plus2 bed site counts, trinuc rates, BisQC
  - Output cols include: `HCG_site_count`, `GCH_site_count`, `chrM_noncpg`, `chr21_noncpg`, `chr21_endo`, `chr21_exo`
- `collect_scnome_qc.py` — aggregates to `scnome_qc_summary.csv`
- `run_qc_and_collect.sh` — single SLURM batch job (loop + collect)

**snmC-seq2** (`snmCseq2/codes/`):
- `snmcseq2_qc_per_cell.py` — parses Bismark SE, BAM summary, CG.6plus2 bed, trinuc rates, BisQC
  - Output cols include: `HCG_site_count`, `chrM_noncpg`, `chr21_noncpg`, `chr21_endo`
  - `--genome hg38|mm10` flag stores genome label
- `collect_snmcseq2_qc.py` — aggregates to `snmcseq2_qc_summary.csv`
- `run_qc_and_collect.sh` — SLURM batch; maps task_id ≤153 → hg38, >153 → mm10

**Nagano** (`nagano/`):
- `schic_qc.py` — fastp JSON + BAM stats → `qc_summary.csv`

**Droplet HiC** (`droplethic/my_project/`):
- `rupture_qc.py` — CB-tagged BAM + pairs → `qc_stats/SRR27586278_hg38.per_cell.tsv`

### Unified comparison

**`summary/compare_qc_all.py`** — reads all method summaries, harmonizes to canonical columns, outputs:
- `all_methods_qc.base_alignment.pdf`
- `all_methods_qc.methylation.pdf`
- `all_methods_qc.hic_contacts.pdf`
- `all_methods_qc.harmonized.csv`
- `all_methods_qc.summary.csv`

Config: `summary/datasets_all.csv`
Runner: `summary/run_compare_qc_all.sh`

**Trinuc pre-computed summaries** (`summary/trinuc/`):
- `scnome.chr21.txt`, `snmCseq2.chr21.txt`, `snmCseq2.chr19.txt`, `smallwood.chr21.txt`, `snmCseq3.chr21.txt`
- Format: `sample\tnoncpg\tendo\texo` (ACT%, ACG%, GCT% respectively)
- For scnome/snmCseq2 files: sample IDs are `{prefix}_1` / `{prefix}_2` (R1/R2 separate)
- `extract_trinuc.py` — generates these from per-cell raw trinuc files
- `compare_qc_all.py` loads these and merges `NonCG_rate` and `GCH_rate` into harmonized data

---

## Key QC Metric Definitions

| Canonical column | Description | Source |
|---|---|---|
| `InputReads` | Total input reads | trim report or Bismark |
| `MapQ30Rate` | % reads with MAPQ≥30 | BAM summary |
| `DupRate` | Duplication rate (%) | BAM summary |
| `FinalMappedReads` | MAPQ30 unique mapped reads | BAM summary |
| `mCG_Rate` | CpG methylation rate (%) | Bismark report |
| `mCH_Rate` | Non-CpG methylation rate (%) | Bismark report |
| `LambdaConvProxy` | Lambda DNA or chrM non-CG meth (%) | BAM or trinuc |
| `NonCG_rate` | ACT% on chr21 (bisulfite conversion proxy) | `summary/trinuc/` |
| `GCH_rate` | GCT% on chr21 (NOMe accessibility) | `summary/trinuc/` |
| `HCG_site_count` | Rows in HCG.6plus2.bed (covered CpG sites) | per-cell bed counting |
| `GCH_site_count` | Rows in GCH.6plus2.bed (covered GCH sites, NOMe) | per-cell bed counting |
| `TotalContacts` | Hi-C contacts | YAP / pairs file |
| `CisLongRatio_Pct` | Cis long-range ratio (%) | YAP / pairs |

**Trinuc context conventions:**
- `noncpg` = ACT% (non-CG, non-GCH — pure conversion proxy)
- `endo` = ACG% (endogenous CpG methylation proxy)
- `exo` = GCT% (GCH, NOMe nucleosome accessibility proxy)
- `HCG count` = number of rows in `HCG.6plus2.bed` (covered CpG sites per cell)
- `GCH count` = number of rows in `GCH.6plus2.bed` (covered GCH sites per cell)

---

## File Naming Conventions

**scNOMe** alignment dir: `04.alignment/`
- Bismark reports: `{prefix}.{prefix}_{1,2}_trimmed_bismark_bt2_SE_report.txt`
- BAM summary: `{prefix}_{1,2}.summary.txt`
- Trinuc: `{prefix}_{1,2}.rmdup.RG.trinuc_methy.{chrM,chr21}.txt`
- HCG bed: `{prefix}_{1,2}.rmdup.RG.cytosine.filtered.sort.HCG.6plus2.bed`
- GCH bed: `{prefix}_{1,2}.rmdup.RG.cytosine.filtered.sort.GCH.6plus2.bed`

**snmC-seq2** alignment dir: `05.align/`
- Bismark reports: `{prefix}_{1,2}.clean_bismark_bt2_SE_report.txt`
- BAM summary: `{prefix}_{1,2}.summary.txt`
- Trinuc: `{prefix}_{1,2}.rmdup.RG.trinuc_methy.{chrM,chr21}.txt`
- CG bed: `{prefix}_{1,2}.rmdup.RG.cpg.filtered.sort.CG.6plus2.bed`

**Smallwood** alignment dir: `05.align_mm10/`
- Trinuc: `{prefix}.rmdup.RG.trinuc_methy.{chrM,chr19}.txt`
- CG bed: `{prefix}.rmdup.RG.cpg.filtered.sort.CG.6plus2.bed`

---

## Recent Changes (last session)

1. **Bisulfite conversion QC clarification:**
   - `HCG_site_count` / `GCH_site_count` = 6plus2 bed row count (covered sites), NOT trinuc context counts
   - Non-CpG methylation rate sourced from `summary/trinuc/` pre-computed files
   
2. **Updated per-cell scripts:**
   - `scnome/codes/scnome_qc_per_cell.py`: bed row counting for HCG + GCH sites; simple ACT/ACG/GCT rate extraction from per-cell trinuc files
   - `snmCseq2/codes/snmcseq2_qc_per_cell.py`: bed row counting for CG sites; ACT/ACG rates from trinuc
   - `smallwood/smallwood_qc_per_cell.py`: simplified trinuc to `_get_trinuc_rate()` helper

3. **Updated `summary/compare_qc_all.py`:**
   - Added `CellID` preservation in harmonized data
   - Added `load_trinuc_data()` + `merge_trinuc_into_hdf()` functions
   - Loads `summary/trinuc/` files; strips `_1`/`_2` and averages R1+R2 for scnome/snmCseq2
   - Added `NonCG_rate` and `GCH_rate` to METHYLATION metric group
   - `harmonize_scnome`: adds `HCG_site_count`, `GCH_site_count`
   - `harmonize_snmcseq2`: adds `HCG_site_count`
   - `harmonize_smallwood`: maps `CpG_TotalSites` → `HCG_site_count`

---

## How to Re-run QC

```bash
# scNOMe
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
sbatch scnome/codes/run_qc_and_collect.sh

# snmC-seq2
sbatch snmCseq2/codes/run_qc_and_collect.sh

# Final comparison plots
cd summary
source /home/jmj7858/.bashrc && conda activate scnomehic
python compare_qc_all.py --config datasets_all.csv --output all_methods_qc
```

---

## Retrieval Command for New Chat

Paste this into a new Cursor chat to restore full context:

```
Read /home/jmj7858/epifluidlab/workspace/scnomehic_paper/benchmark/PROJECT_CONTEXT.md and use it as the full project context. I am continuing work on the scNOME-HiC benchmark project. The key files are in /home/jmj7858/epifluidlab/workspace/scnomehic_paper/benchmark/. Please confirm you have read and understood the context before proceeding.
```
