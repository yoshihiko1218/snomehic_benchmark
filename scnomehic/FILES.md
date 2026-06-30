# scnomehic — file/folder guide

scNOMe-HiC — **this is the published method**. GM12878 cell line, hg38, 188 cells
(`acc_list.txt`, `scNH_GM_4plex_*`). Produces three modalities: CpG methylation
(HCG), GpC accessibility (GCH), and Hi-C contacts.

## Real run path (the benchmark)
Two alignment tracks; the **bhmem track is the method**, the **YAP track is the
local comparison aligner**:
- **YAP / m3C (local, canonical for QC here):** TrimGalore → `yap start-from-cell-fastq`
  (`mode=m3c`, bismark + **bowtie2**) → allcools methylation + m3C `generate-contacts`.
  Output → `alignment/` (255G). **`alignment/stats/MappingSummary.csv.gz` is read by
  `summary/qc.ipynb`** (Hi-C contacts + mapping metrics).
- **bhmem (the method aligner) — output is EXTERNAL, not in this folder:**
  `bisulfitehic` BWA-MEM bisulfite aligner. The valid `.calmd.bam` / trinuc / HCG+GCH
  methylation live at
  `/projects/b1198/epifluidlab/yoshii/sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake/`
  and `.../07.bistools_snakemake/methylation/` — these feed the scnomehic HCG/GCH
  columns in `summary/` (`gch_hcg_counts/scnomehic*.gchhcg.txt`).
  NOTE: there is **no valid local bhmem BAM dir** here — an earlier local bhmem run
  failed (munged `_1`/`_2` read IDs in `fastq/` broke BWA pairing; see JOBS.md), so
  bhmem was rerun on the external b1198 copy.

### Keep
| Item | Size | Role |
|---|---|---|
| `alignment/` | 255G | YAP/bowtie2 — `stats/MappingSummary.csv.gz` consumed by qc.ipynb |
| `fastq/` | 113G | raw input (read IDs munged; fine for YAP, not for bhmem) |
| `codes/`, `logs/`, `acc_list.txt`, `JOBS.md`, `SESSION_NOTE_2026-05-28_1.md` | small | pipeline + records |

## Dead-end experiment (safe to delete — unreferenced by any QC pipeline)
Verified 2026-06-29: not read by `summary/` or `qc.ipynb`.
- **`alignment_bowtie1/` (341G)** — bowtie1 vs bowtie2 parameter experiment
  (different clip lengths, separate bismark index). Superseded by `alignment/`
  (bowtie2), which is the canonical local arm.
