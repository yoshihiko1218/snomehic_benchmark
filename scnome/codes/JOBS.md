# JOBS — scnome

## GM12878 rerun (2026-06-09) — both-ends trim fix, NOT YET SUBMITTED
Reason: GM12878 must clip 6 bp from BOTH ends (Pott 2017); previous run used 5'-only.
acc_list.txt trimmed to 34 cells (controls SRR3729654-3729660 removed; backup
`acc_list_with_controls.txt`). All step arrays scoped to `--array=1-12` (GM cells).
Deletions DONE (2026-06-09): stale GM outputs + control files cleared.
Note: raw reads are gzipped `01.fastq/<p>_{1,2}.fq.gz`; 01.trim.sh now reads
`.fq.gz` (was `.fastq`), and scnome_qc_per_cell.py accepts either trim-report name.

Submit the whole chain with ONE command:
- `bash codes/submit_gm.sh`
  -> 01.trim -> 02.alignment -> 03.methy_extract -> run_qc (each afterok the prev).

SUBMITTED 2026-06-09 (array 1-12 = GM12878):
- 01.trim.sh          : 4257035
- 02.align_dedup.sh     : 4257036  (afterok 4257035)
- 03.methy_extract.sh : 4257037  (afterok 4257036)
- 04.qc_per_cell.sh           : 4257038  (afterok 4257037)
=> FAILED. Root cause: 01.trim.sh `cd` pointed at nonexistent b1198 path, so
   acc_list.txt was not found, prefix empty, trim_galore ran on missing files
   yet exited 0; 02.alignment then ran with no trimmed input and FAILED; 03/qc
   DependencyNeverSatisfied (cancelled 4257037/4257038). 01.trim's 25s runtime
   was the tell. Fixed cd (b1198->b1042) + added fail-fast guards.

RESUBMITTED 2026-06-09 (array 1-12 = GM12878):
- 01.trim.sh          : 4258488
- 02.align_dedup.sh     : 4258489  (afterok 4258488)
- 03.methy_extract.sh : 4258490  (afterok 4258489)
- 04.qc_per_cell.sh           : 4258491  (afterok 4258490)
Logs: logs/01.qc_trim/, logs/02.bisqc/, logs/03.methy_extract/, logs/03.qc/
Monitor: squeue -u jmj7858

STATUS 2026-06-09 ~18:21:
- 01.trim (4258488): all 12 COMPLETED (3-9 min each). 24 trimmed .fq.gz present.
  GM both-ends clip CONFIRMED in report (6 bp 5' AND 6 bp 3').
- 02.alignment (4258489): all 12 RUNNING, bismark progressing, no errors.
- 03.methy (4258490) / run_qc (4258491): PENDING (dependency).
- Note: 28-byte *.rmdup.bam at 17:39 are stale leftovers from the failed first
  chain; current run overwrites them after bismark. (02.align_dedup.sh has no
  fail-fast; harmless here since bismark is succeeding.)

STATUS 2026-06-09 ~19:50:
- 02.alignment (4258489): 2/12 COMPLETED (tasks 3,5 = SRR3729644/646), 10 RUNNING
  (~1h40m). Completed cells produced REAL BAMs (256-350 MB .rmdup.bam +
  .rmdup.RG.bam + .bai + summary.txt); stale 28-byte files overwritten. No errors.
  Non-directional library confirmed (CT/CT, CT/GA, GA/CT, GA/GA all present).
- 03.methy (4258490) starts only after ALL 12 alignment tasks finish (afterok).

STATUS 2026-06-09 ~21:59:
- 02.alignment (4258489): ALL 12 COMPLETED (exit 0). 12/12 GM cells have all 4
  real BAMs (_1/_2 x rmdup/rmdup.RG, >1MB). Stale 28-byte files all overwritten.
- 03.methy_extract (4258490): RUNNING all 12 (~11 min). bismark_methylation_extractor
  producing .bismark.cov.gz/bedGraph/M-bias/splitting reports. No errors yet.
  Next: coverage2cytosine --nome-seq -> NOMe.{CpG,GpC} reports.
- run_qc (4258491): PENDING (afterok 03).

FINAL STATUS 2026-06-10 ~00:50 — GM12878 RERUN COMPLETE (all 4 jobs exit 0):
- 01.trim     (4258488): 12/12 COMPLETED. GM both-ends clip applied (6bp 5'+3').
- 02.alignment(4258489): 12/12 COMPLETED. All cells: real _1/_2 rmdup(.RG).bam.
- 03.methy    (4258490): 12/12 COMPLETED. 24 NOMe.CpG.cov.gz + 24 NOMe.GpC.cov.gz.
- run_qc      (4258491): 12/12 COMPLETED (~26s each). 12 qc_stats/*.qc_stats.csv.
Verified: 24 trimmed.fq.gz, 12/12 cells x4 real BAMs, 24+24 NOMe cov.gz, 12 QC csv.
No errors (only benign 'could not extract chromosomal sequence ..._random').
NOTE: QC columns sourced from the OLD Bis-tools route (HCG/GCH_site_count,
chrM/chr21 trinuc rates) are blank for these cells - that route is not part of
the Bismark pipeline. Core Trim/Bismark/BAM/methyl metrics are populated.
NEXT (optional): `sbatch codes/04.qc_collect.sh` to rebuild
scnome_qc_summary.csv across ALL 34 cells (GM rerun + existing K562).

## K562 per-cell pipeline (merge 2 runs/cell) — SUBMITTED 2026-06-10
K562 cells = 11 (consecutive SRR pairs, acc_list_k562_cells.tsv). Each cell was
sequenced as 2 runs; merge both runs' bismark BAMs then markdup ONCE (removes
within+cross-run dups), then NOMe methylation per cell. Reuses existing per-run
bismark BAMs (skips trim+align).
- 02k.merge_dedup_k562.sh   : 4288785  (array 1-11)
- 03k.methy_extract_k562.sh : 4288786  (afterok 4288785, array 1-11)
Submit driver: `bash codes/submit_k562.sh`. Logs: logs/04.k562_merge/, logs/04.k562_methy/.

## NOMe QC (trinuc + detected HCG/GCH sites) — Bismark-native, DONE/SUBMITTED 2026-06-10
codes/nome_qc_sites_trinuc.py + codes/04.qc_nome_sites.sh. Computes per cell:
HCG/GCH detected sites (NOMe cov rows) + chrM/chr21 ACT/ACG/GCT meth rates
(from deduped BAM + reference). Validated on GM SRR3729642: chr21 ACT=6.5%
(conversion), ACG=38%, GCT=26%; chrM ~75% = known mito bisulfite-resistance.
- K562 merge+dedup (4288785): 11/11 COMPLETED (~7 min). 11/11 real merged BAMs.
- K562 methy (4288786): RUNNING (array 1-11).
- GM NOMe QC (4289055): array 1-12, LISTFILE=acc_list.txt -> qc_stats/<SRR>.nome_qc.tsv.
- K562 NOMe QC (4289056): array 1-11, afterok 4288786 -> qc_stats/K562_NN.nome_qc.tsv.
Submit form: sbatch --array=1-N --export=ALL,LISTFILE=<list> codes/04.qc_nome_sites.sh.

## BisSNP trinuc QC (ORIGINAL method) — SUBMITTED 2026-06-10
Decision: use the original BisSNP for trinuc (my read-level proxy differed: it
lacked the -minPatConv 0.8 read conversion filter, used context CLASSES not
literal trinucleotides, and no dbSNP masking -> read several-fold higher).
codes/04.qc_bissnp_trinuc.sh: BisulfiteGenotyper on chrM+chr21 per cell/mate ->
04.alignment/<cell>_<m>.rmdup.RG.trinuc_methy.{chrM,chr21}.txt (consumed by
scnome_qc_per_cell.py). CRITICAL: needs Java 8 (module java/jdk1.8.0_191);
GATK-3.8 walker discovery fails on the conda env's Java 21. Validated on
SRR3729642 chrM (exit 0, 16 N-C-N contexts).
- GM trinuc   (4290691): array 1-12, LISTFILE=acc_list.txt
- K562 trinuc (4290692): array 1-11, LISTFILE=acc_list_k562_cellids.txt
After these finish: run_qc / run_qc_and_collect populates trinuc columns.
(Note: nome_qc HCG/GCH = covered-site counts from NOMe cov, NOT the BisSNP
6plus2-filtered counts; different definition.)


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
