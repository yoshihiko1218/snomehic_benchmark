# Bisulfite-conversion figures + snmCAT base-alignment metrics

Added 2026-06-30. Adds the missing **bisulfite-conversion** distribution (metric 1)
across the 6 methylation methods, and folds **snmCAT** into the readcount/mapq30
panels (metrics 2 & 3). All wired into `qc.ipynb` and regenerated end-to-end.

## Conversion proxy
`noncpg = ACT%`, `endo = ACG%`, `exo = GCT%` (trinucleotide methylation rates).
`noncpg` (non-CpG, non-GpC) = bisulfite non-conversion proxy. Reported for two
contexts: **chrM** and an **autosome** (chr21 human / chr19 mouse).

## Scripts (run in this order; conda env `scnomehic`, from benchmark root)
1. `extract_trinuc.py --folder snmCseq3/04.bhmem_bam --suffix .calmd.trinuc_methy.chrM.txt --output summary/trinuc/snmCseq3.chrM.txt`
   — snmCseq3 chrM table (chr21 already existed).
2. `python summary/extract_trinuc_snmCAT.py`
   — snmCAT had no BisSNP/bhmem trinuc files; computes ACT/ACG/GCT % per cell from
   the YAP **allc** files (tabix-indexed) for chrM + chr21. 4-mer allc context →
   trinuc = first 3 chars; rate = Σmc/Σcov (position-level). Writes
   `summary/trinuc/snmCAT.{chrM,chr21}.txt`.
3. `python summary/collect_conversion.py`
   — assembles `summary/conversion_percell.csv` (554 cells; cols dataset, cell,
   noncpg_chrM, noncpg_auto, auto_chrom). **Per-cell, never per-mate.** Cell sets
   match the other panels: scnome 23, smallwood 51, snmCseq2 96 (mm10), snmCseq3 98,
   scnomehic 187 (passed), snmCAT 99. Sources: per-cell qc_summary chrM/chr-auto for
   scnome/smallwood; snmCseq2 chrM from qc_summary + chr19 from per-mate trinuc
   averaged; trinuc tables for snmCseq3/snmCAT; external `gm.{chrM,chr21}` for scnomehic.
4. `python summary/plot_conversion.py` (standalone) OR the appended qc.ipynb cells
   — write `figures/conversion_{chrM,chr21chr19}_violin.{pdf,png}`.

## qc.ipynb changes
- New snmCAT-build cell (after `df_all_gm`): `mapping_efficiency_pct =
  mean(R1,R2 MappingRate)`, `unique_best_hit_n = mean(R1,R2 UniqueMappedReads)`,
  from `snmCAT/mapping_brain/stats/MappingSummary.csv.gz`, restricted to 99 cells.
- snmCAT added to the readcount (metric 2) + mapq30 (metric 3) violins.
- Two appended cells plot the conversion figures from `conversion_percell.csv`.

## Caveat
snmCAT's MappingSummary has **no MapQ30 columns** (unlike snmCseq3), so its
"% uniquely-mapped" is the unique mapping rate (mean R1/R2 `MappingRate`), NOT
MapQ30-filtered — it therefore reads slightly higher than the MapQ30-based methods.
snmCAT conversion is position-level (allc Σmc/Σcov); the others are read-level ACT%.
Same ACT/ACG/GCT proxy, slightly different weighting.
