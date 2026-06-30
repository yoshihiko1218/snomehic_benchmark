# snmCseq3 — file/folder guide

snm3C-seq benchmark (mm10, ~100 cells). This folder doubled as the arena for a
deep **bhmem-vs-YAP alignment methodology investigation** (NM correction, XG/XR
tags, MAPQ), so it accumulated many experiment dirs. This guide separates the
**canonical benchmark run** (consumed by `summary/qc.ipynb`) from **dead-end
experiments**.

## Canonical run path (the benchmark)
`fastq/` (raw) → trim → `03.trimmed_fastq/` → two alignments, **both consumed by
`summary/qc.ipynb`**:
- **YAP** → `alignment/stats/MappingSummary.csv.gz` — Hi-C contacts + mapping
  metrics (qc.ipynb lines ~3061/3105; also `summary/datasets_all.csv`).
- **bhmem** → `04.bhmem_bam/` — `*.summary.txt` alignment metrics (qc.ipynb
  ~2543), `*.calmd.trinuc_methy.*` (→ `summary/trinuc/`), `*.CG.6plus2.bed`
  (→ HCG counts in `summary/count_*`).

### Keep
| Item | Role |
|---|---|
| `fastq/` | raw input |
| `03.trimmed_fastq/` | trimmed input (regenerable from fastq) |
| `alignment/` | YAP mapping — consumed by qc.ipynb |
| `04.bhmem_bam/` | bhmem BAMs/summaries — consumed by qc.ipynb |
| `02.fastqc_out/` | FastQC reports |
| `codes/` | pipeline + `codes/comparison/` investigation scripts + method `.md` docs |
| `logs/`, `acc_list*.txt`, `SESSION_NOTE_2026-04-04_1.md` | run records |
| `final_comparison/`, `comprehensive_comparison/`, `corrected_nm_comparison_v2/` | small definitive investigation results (records) |
| `bhmem_original_vs_xg_comparison.tsv`, `bhmem_old_vs_new_comparison.tsv` | tiny conclusion tables |

## Cell counts
| Stage | Count |
|---|---|
| Listed in `acc_list.txt` (full SRA pool) | 1379 |
| Actually aligned (`04.bhmem_bam/*.calmd.bam`; YAP `MappingSummary` rows) | 100 |
| **Used in benchmark QC figures** (`.summary.txt` + trinuc + valid `trans_ratio`) | **98** |

`acc_list_paired.txt` (61) and the `*_subset` dirs were investigation batches, not
the QC set.

## QC-metric provenance (what feeds each benchmark distribution)
All metrics derive from small text outputs, NOT the heavy `.calmd.bam`. Both KEEP
dirs (`04.bhmem_bam/`, `alignment/`) are required by `summary/qc.ipynb`.

Verified against `qc.ipynb` execution order: cell 29 first builds `snmCseq3`
from bhmem `*.summary.txt`, but cells 32/33 **overwrite** `snmCseq3` from the YAP
`MappingSummary.csv.gz` before any plot runs. So the **final figures source
alignment + Hi-C contacts from YAP**, and **conversion + HCG-loci from bhmem**:

| Benchmark distribution | Tool (final figure) | File read |
|---|---|---|
| Bisulfite conversion (chrM, chr19) | **bhmem** | `04.bhmem_bam/*.calmd.trinuc_methy.{chrM,chr19}.txt` → `summary/trinuc/snmCseq3.chr21.txt` (cell 28) |
| Count uniquely-mapped fragments | **YAP** | `alignment/stats/MappingSummary.csv.gz` (`R1/R2UniqueMappedReads` → `mapped_read_n`, cell 33/48) |
| % uniquely-mapped (MAPQ30 rate) | **YAP** | `MappingSummary` (`R1/R2MappingRateMapQ30` → `mapping_rate_mapq30`, cell 47) |
| per-cell cis-contacts | **YAP** | `MappingSummary` (`CisShortContact + CisLongContact`, cell 49) |
| Trans/Cis ratio | **YAP** | `MappingSummary` (`TransRatio`, cell 49) |
| per-cell cis-contacts >1kb | **YAP** | `MappingSummary` (`CisLongContact`, cell 49) |
| (HCG loci — separate figure) | **bhmem** | `04.bhmem_bam/*.CG.6plus2.bed` → `summary/gch_hcg_counts/` |

Net: **both dirs are required** — YAP `alignment/` for align+contacts (metrics
2–6), bhmem `04.bhmem_bam/` for conversion (metric 1) + HCG loci. The bhmem
`*.summary.txt` contact columns are computed in cell 29 but unused (overwritten),
so bhmem is NOT the contact source despite using `mh_reads_summary.v2.py`.

## Dead-end experiments (safe to delete — unreferenced by any QC pipeline)
Verified 2026-06-29: none of the below are read by `summary/` scripts or
`qc.ipynb`; only `tools/commit_all.sh` git-adds `alignment_bowtie1`/`_mapq0`
(those lines removed during cleanup).

`alignment_mapq0/`, `alignment_bowtie1/`, `04.bhmem_bam_subset/`,
`04.bhmem_bam_xg/`, `04.bhmem_bam_original/`, `04.bhmem_bam_oldjar/`,
`04.bhmem_bam_buf100k/`, `04.bhmem_bam_noxg/`, `test_fastq/`,
`test_fastq_renamed/`, `03.trimmed_fastq_subset/`, `mapq_comparison/`,
`mapq_comparison_mapq0/`, `corrected_nm_comparison/` (round-1, superseded by
`_v2`), and loose root scratch BAMs `one.*`, `pdf_subset.*`.

Investigation conclusions: original & XG JARs produce identical alignments; the
Feb R2 low-MAPQ issue was input-FASTQ-dependent, not a code bug; corrected-NM
method validated to 99.5% (see `codes/comparison/BISULFITE_NM_CORRECTION_METHOD.md`
and `SESSION_NOTE_2026-04-04_1.md`).
