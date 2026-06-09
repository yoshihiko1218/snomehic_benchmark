# JOBS — scnome

## GM12878 rerun (2026-06-09) — both-ends trim fix, NOT YET SUBMITTED
Reason: GM12878 must clip 6 bp from BOTH ends (Pott 2017); previous run used 5'-only.
acc_list.txt trimmed to 34 cells (controls SRR3729654-3729660 removed; backup
`acc_list_with_controls.txt`). All step arrays scoped to `--array=1-12` (GM cells).
Deletions DONE (2026-06-09): stale GM outputs + control files cleared.
Note: raw reads are gzipped `01.fastq/<p>_{1,2}.fq.gz`; 01.trim.sh now reads
`.fq.gz` (was `.fastq`), and scnome_qc_per_cell.py accepts either trim-report name.

Submit the whole chain with ONE command:
- `bash codes/submit_gm_rerun.sh`
  -> 01.trim -> 02.alignment -> 03.methy_extract -> run_qc (each afterok the prev).

SUBMITTED 2026-06-09 (array 1-12 = GM12878):
- 01.trim.sh          : 4257035
- 02.alignment.sh     : 4257036  (afterok 4257035)
- 03.methy_extract.sh : 4257037  (afterok 4257036)
- run_qc.sh           : 4257038  (afterok 4257037)
Logs: logs/01.qc_trim/, logs/02.bisqc/, logs/03.methy_extract/, logs/03.qc/
Monitor: squeue -u jmj7858


## 03.methy_extract — Bismark NOMe methylation extraction (Pott 2017 protocol)
- **Job name:** methyext
- **Job ID:** 4128017  (array 1-41, one task per cell in `acc_list.txt`)
- **Submitted:** 2026-06-08
- **Command:** `sbatch codes/03.methy_extract.sh`
- **Logs:** `scnome/logs/03.methy_extract/methyext.<arrayid>.txt` (+ `.err`)
- **What it does:** for each cell's `_1`/`_2` `rmdup.bam`:
  - Step 1 `bismark_methylation_extractor -s --ignore 6 --bedGraph --CX` →
    `05.methy/<cell>_<mate>.rmdup.bismark.cov.gz` (ALREADY PRESENT for most
    cells, so this step is skipped on resume).
  - Step 2 `coverage2cytosine --nome-seq` →
    `05.methy/<cell>_<mate>.NOMe.CpG.cov.gz` (ACG/TCG, CpG) and
    `05.methy/<cell>_<mate>.NOMe.GpC.cov.gz` (GCA/GCC/GCT, GpC); ambiguous GCG
    dropped. Per-locus, per-strand, no SNP filtering — matches the YAP/allcools
    convention (NOT BisSNP).
- **Notes:** coverage2cytosine OOMs on the login node (loads whole genome);
  must run via SLURM (64 G). Resume-safe: skips a mate when both NOMe reports
  exist.
- **Next:** after completion, count detected loci with
  `summary/count_bismark_nome_loci.py` (rows in each `.NOMe.{CpG,GpC}.cov.gz`).

### Status checks
- `squeue -u jmj7858 | grep methyext`
- On failure inspect `logs/03.methy_extract/methyext.<id>.err`.
