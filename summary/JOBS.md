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
