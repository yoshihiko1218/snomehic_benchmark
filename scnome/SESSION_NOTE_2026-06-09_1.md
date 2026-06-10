# Session Note — 2026-06-09 (seq 1)

Folder: benchmark/scnome

## Goal
Verify scNOMe-seq trimming/alignment code vs Pott 2017 methods, fix GM12878
both-ends clipping, scope a GM12878 rerun, and remove control samples.

## Findings
- Trim/align params matched the paper EXCEPT the GM12878-specific rule: "6 bp
  clipped from either end" was not implemented (code only had `--clip_R1 6`).
- This data: SRR3729642-3729653 = GM12878, SRR3729654-3729660 = controls
  (spike-in, exclude), SRR3729661-3729682 = K562.
- GM outputs already existed but were built with old 5'-only trimming → STALE.
  Trimmed `.fq.gz` already deleted (only reports remain).
- `02.alignment.sh` had bismark/dedup/RG steps all commented out (only
  bam_summary active) → could not re-align as-is.

## Changes made
- `codes/01.trim.sh`: classify by SRR; GM → `--clip_R1 6 --three_prime_clip_R1 6`,
  K562 → `--clip_R1 6`, controls → skip. Array 1-41 → 1-12.
- `codes/02.alignment.sh`: control-skip guard; ENABLED full chain bismark →
  sort|markdup → addreplacerg(+index) → bam_summary. Array 1-12. markdup kept
  WITHOUT `-r` (same dedup as all cell types; only trimming differs for GM).
- `codes/03.methy_extract.sh`, `codes/run_qc.sh`: array 1-12 + clarifying comments.
- `acc_list.txt`: removed 7 controls → 34 cells (GM 1-12, K562 13-34). Backup:
  `acc_list_with_controls.txt` (41 lines).
- New `codes/clear_gm_stale.sh` (dry-run/`--yes`): deletes stale GM outputs (842 files).
- New `codes/clear_control_files.sh` (dry-run/`--yes`): deletes control files (258 files).
- New `codes/FILES.md`; updated `codes/JOBS.md` with rerun plan.

## Verification
- `bash -n` passes on all scripts.
- Classification logic checked across all 41 original accessions (12/7/22).
- New acc_list.txt = 34 cells, no controls.
- Dry-runs: GM stale = 842 files, controls = 258 files.

## Not done / next (user-driven)
- Deletions NOT executed (no delete permission). User runs:
  `bash codes/clear_gm_stale.sh --yes` and `bash codes/clear_control_files.sh --yes`.
- GM rerun jobs NOT submitted yet (see codes/JOBS.md submit order).

## GM rerun execution
- Submitted chain 4257035-4257038 -> FAILED. Cause: 01.trim.sh `cd` used wrong
  project path (b1198 vs b1042); acc_list.txt not found, empty prefix, trim ran
  on missing files but exited 0; 02.alignment FAILED on no input; 03/qc never ran.
- Also fixed earlier: 01.trim referenced `.fastq` inputs but raw reads are
  `.fq.gz`; QC parser now accepts either trim-report name.
- Fixed cd path + added fail-fast guards (bad cd / empty prefix / missing input /
  trim failure / missing output). Cancelled 4257037,4257038.
- Resubmitted chain 4258488(trim)->4258489(align)->4258490(methy)->4258491(qc).
  01.trim now RUNNING all 12 tasks (no longer instant-failing). Monitoring.

## GM rerun COMPLETE (2026-06-10 ~00:50)
- Chain 4258488->4258489->4258490->4258491 all COMPLETED (exit 0).
- Mid-run fixes applied: (1) 01.trim cd b1198->b1042 + fail-fast guards;
  (2) 01.trim input .fastq->.fq.gz (raw reads are gzipped); (3) QC parser accepts
  either .fq.gz/.fastq trim-report name.
- Verified outputs (GM SRR3729642-3729653): 24 trimmed.fq.gz; 12/12 cells x4
  real rmdup(.RG).bam; 24 NOMe.CpG.cov.gz + 24 NOMe.GpC.cov.gz; 12 qc_stats csv.
- GM both-ends clip (6bp 5'+3') confirmed in trimming reports. No real errors.
- Known: QC columns from the old Bis-tools route are blank (route not in Bismark
  pipeline). Optional next: run_qc_and_collect.sh to rebuild summary over 34 cells.
- Monitoring stopped.

## ALL COMPLETE (2026-06-10)
- GM rerun (both-ends clip): 12 cells through trim->align->dedup->methy->QC.
- K562: 11 merged cells (consecutive SRR pairs) via merge-2-runs -> dedup once ->
  methy. K562 methy (4288786) + NOMe QC (4289056) both 11/11 COMPLETED.
- NOMe cov: GM 24+24, K562 22+22 (CpG+GpC). nome_qc.tsv: GM 12 + K562 11.
- BisSNP trinuc (original method; needs Java 8 module jdk1.8.0_191): 46 chrM +
  46 chr21 trinuc_methy files. Validated; chr21 conv ~0.3% (GM).
- Compared to gm_sc_new (scNOMe-HiC, 188 cells): conversion & CpG comparable;
  GpC accessibility much higher in gm_sc_new (~65% vs ~18%) = stronger enzyme/
  accessibility signal in the in-house scNOMe-HiC protocol.
- Renamed codes/ by execution stage (01/02/02k/03/03k/04.qc_*, submit_*, util_*).
- Optional remaining: 04.qc_collect.sh to rebuild scnome_qc_summary.csv.
- Monitoring stopped.

## Commits
- cfd2cc1 — per-cell-type clipping + exclude controls (earlier).
- (this session) — GM rerun scoping, alignment chain enable, acc_list trim, helpers, docs.
