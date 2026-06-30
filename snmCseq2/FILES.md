# snmCseq2 — file/folder guide

snmC-seq2 benchmark. **Mixed species: 249 cells = 153 hg38 + 96 mm10** (each cell's
`_1`/`_2` mates share one genome; `snmCseq2_genome_map.tsv` assigns genome).
Not NOMe → CpG methylation only, no GCH.

## Real run path — TWO arms were run
1. **Bismark-SE per-mate (canonical; feeds the benchmark figures):**
   `01.fastq` → cutadapt (`03.trimmed_fastq`) → `bismark -bowtie2` SE (R1 `-pbat`,
   R2 `-non_directional`) → `05.align/*.clean_bismark_bt2.bam` → `deduplicate_bismark`
   + `bismark_methylation_extractor` → `06.methy/*.cov.gz` → `snmcseq2_qc_per_cell.py`
   → **`snmcseq2_qc_summary.csv`** (249 rows). HCG loci destranded from
   `06.methy/*.deduplicated.bismark.cov.gz` (`summary/count_hcg_destranded.py` ds=snmCseq2).
2. **YAP / cemba_data `mc` (the user's parallel YAP run on the same cells):**
   `yap_mapping_hg38/` (153 hg38 cells) + `yap_mapping_mm10/` (96 mm10 cells), each
   → `stats/MappingSummary.csv.gz` + `AllcPaths.tsv`. Used as a Bismark-vs-YAP HCG
   sanity arm (`count_hcg_destranded.py` ds=`snmCseq2_yap_mm10`, reads
   `yap_mapping_mm10/stats/AllcPaths.tsv`). KEEP — deliberately produced.

### What the summary page consumes (provenance)
| Benchmark distribution | Tool | File |
|---|---|---|
| Bisulfite conversion (chrM, chr21/chr19) | Bismark (trinuc) | `summary/trinuc/snmCseq2.chr{21,19}.txt` (from `05.align/*.rmdup.RG.trinuc_methy.*`) |
| Count / % uniquely-mapped | Bismark | `snmcseq2_qc_summary.csv` (`unique_best_hit_n`, `mapping_efficiency_pct`) |
| HCG loci (hcg_loci figure) | Bismark cov | `06.methy/*.deduplicated.bismark.cov.gz` → `summary/gch_hcg_counts/` |
| (HCG sanity cross-check) | YAP allc | `yap_mapping_mm10/stats/AllcPaths.tsv` |

snmC-seq2 is **not** a Hi-C method → not in the cis/trans contact figures.

## Cell counts
| Stage | Count |
|---|---|
| Cells (acc_list / genome_map) | 249 (153 hg38 + 96 mm10) |
| QC summary rows | 249 |

## Cleanup candidates
### Big intermediate (safe to delete — ~612G)
`06.methy/{CpG,CHG,CHH}_context_*.txt` — per-read context dumps emitted by
`bismark_methylation_extractor`. Verified 2026-06-29 unreferenced; the distilled
`.cov.gz` (498, consumed for HCG), `.bedGraph.gz`, splitting/M-bias/dedup reports,
`snmcseq2_qc_summary.csv`, and `summary/trinuc/snmCseq2.*` all already exist.
(CHH 447G + CHG 136G + CpG 30G.)

### Small dead-ends
- `.old_single_genome_attempt/` (1.7M) — abandoned mm10-only YAP run (job 4327205,
  cancelled when mixed-species discovered). 0 refs.
- `04.fastqc_out_2/` (313M) — post-trim FastQC reports; regenerable, not consumed (optional).
- `fastq_yap_hg38/`, `fastq_yap_mm10/` (small) — symlink staging dirs for the YAP run.

### Keep (canonical / consumed / deliberately produced)
`01.fastq`, `03.trimmed_fastq`, `05.align`, `06.methy` (minus context txt),
`yap_mapping_hg38`, `yap_mapping_mm10`, `snmcseq2_qc_summary.csv`,
`snmCseq2_genome_map.tsv`, `methy_samples.txt`, `codes/`, `logs/`, `qc_stats/`,
`06.summary/`.
Second-tier (regenerable, kept by default): `05.align` BAMs — trinuc/beds already
extracted; legacy BisSNP `.rmdup.RG.cpg.*.vcf` artifacts also live here.
