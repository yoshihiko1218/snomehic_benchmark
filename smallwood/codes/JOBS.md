# JOBS — smallwood

## 05.hcg_nome — HCG track for cross-method benchmark vs scNOMe
- **Job name:** smwhcg
- **Job ID:** 4299883 (array 1-32, submitted 2026-06-10 ~12:18). Each task ~27 min.
- **Script:** `codes/05.hcg_nome.sh` (array 1-32, reads `acc_list_esc.txt` = 32 ESC cells)
- **Submit:** `sbatch codes/05.hcg_nome.sh`
- **Status check:** `squeue -u jmj7858 -n smwhcg` ; done = `grep -l "DONE prefix" logs/05.hcg_nome/*.txt | wc -l`
- **Logs:** `smallwood/logs/05.hcg_nome/smwhcg.<arrayid>.txt` (+ `.err`)
- **What it does:** `coverage2cytosine --nome-seq --genome_folder mm10_bismark`
  on each `06.methy/<cell>.dedup.bismark.cov.gz` (all-CpG) ->
  `06.methy/hcg/<cell>.NOMe.CpG.cov.gz` = HCG track. Bismark 0.24.2, identical
  tool/version to scNOMe (`scnome/codes/03.methy_extract.sh`), so HCG def is
  identical by construction. Confirmed via `codes/_test_hcg_one.sh`: NOMe-seq
  reports **ACG and TCG context only** (drops GCG and CCG); runs fine on the
  CpG-only cov (no --CX needed).
- **Count step (after array):** `python codes/count_hcg_sites.py --meta
  metadata_esc.tsv --methy 06.methy --out 06.methy/hcg/hcg_site_counts.tsv`
  -> per-cell allCpG vs HCG site counts (HCG = rows of NOMe.CpG.cov.gz, same
  definition as scNOMe `nome_qc_sites_trinuc.py::count_cov_sites`).
- **Caveat (M-bias):** Smallwood 5'-clipped 9 bp at trim; scNOMe used clip 6 +
  extractor `--ignore 6`. Minor stage difference, both remove priming bias.

## 04.methy_extract — Bismark CpG methylation extraction (Smallwood 2014 scBS-seq)
- **Job name:** smwmethy
- **Job ID:** 4128321  (array 1-51, one task per cell in `acc_list.txt`)
- **Submitted:** 2026-06-08
- **Command:** `sbatch codes/04.methy_extract.sh`
- **Logs:** `smallwood/logs/04.methy_extract/smwmethy.<arrayid>.txt` (+ `.err`)
- **What it does:** scBS-seq is NOT NOMe -> CpG only (no GpC/GCH).
  - Step 0: exclude marked duplicates (`samtools view -F 1024`; *.rmdup.bam are
    markdup'd ~21% but not removed) — protocol: "calls extracted after
    duplicates excluded".
  - Step 1: `bismark_methylation_extractor -s --comprehensive --bedGraph
    --genome_folder mm10_bismark` -> `06.methy/<cell>.dedup.bismark.cov.gz`
    (one row per covered CpG; per-locus, per-strand, no SNP filtering — YAP/
    Bismark convention, NOT BisSNP). Detected CpG (HCG) loci = #rows.
  - Temp `<cell>.dedup.bam` removed after the cov is written.
- **Validation:** SRR1248472 (shallow, 32,594 reads) -> 9,791 CpG loci vs its
  BisSNP CG.6plus2 = 1,369 (7.2x more; consistent with per-locus vs BisSNP).
- **Next:** count loci into `summary/gch_hcg_counts/` after completion.

### Status checks
- `squeue -u jmj7858 | grep smwmethy`
- On failure inspect `logs/04.methy_extract/smwmethy.<id>.err`.
