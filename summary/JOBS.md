# JOBS — summary/

## scnome per-cell destranded HCG merge
- **Job ID:** 4893170  (`scnome_hcg_pc`)
- **Submit:** `sbatch summary/submit_scnome_hcg_percell.sh` (from summary/)
- **Account/Partition:** b1042 / genomics, 48G, 8 cpus
- **Script:** `summary/scnome_hcg_percell.py` — unions the two R1/R2 mate covs per
  cell (23 cells) before destranded HCG counting (can't sum per-mate counts).
- **Output:** `summary/gch_hcg_counts/scnome.hcg_percell_destranded.txt`
- **Log:** `summary/logs/hcg_genome/scnome_hcg_percell.txt`
- **Status:** submitted 2026-06-18 (PENDING).

## scnome per-cell HCG+GCH (combined) -- CANONICAL
- **Job ID:** 4893237 (`scnome_loci_pc`), COMPLETED 4:11.
- **Script:** `summary/scnome_loci_percell.py` (+ submit_scnome_loci_percell.sh)
- **Output:** `summary/gch_hcg_counts/scnome.loci_percell.txt` (23 cells)
  median HCG 1,590,594 | GCH 8,676,854 (per-cell union of R1/R2, destranded).
- make_dataset_summaries.py now reads scnome HCG+GCH from this file.
- Supersedes the HCG-only job 4893170 / scnome.hcg_percell_destranded.txt.
