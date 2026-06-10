# SESSION NOTE 2026-06-10_1 — smallwood

## Goal
Decide which Smallwood 2014 scBS-seq cells to benchmark, audit the pipeline
against the paper/GEO methods, and add an HCG track comparable to the scNOMe
benchmark.

## Decisions
- **Benchmark cells = ESCs only:** 12 x 2i (SRR1248457–1248468) + 20 x Ser
  (SRR1248477–1248496) = 32 cells. Excluded CT/oocyte/other (SRR1248444–455,
  470–476) and all bulk/"deeper" samples. Ranges confirmed by user.
- **Methylation calling = Bismark** (GEO `*cov.txt`/`*CpG.txt` formats are exact
  Bismark methylation-extractor / coverage2cytosine outputs). Paper reports
  all-CpG; no separate caller.
- **HCG for benchmark:** generate via `coverage2cytosine --nome-seq` (same tool
  Bismark 0.24.2 as scNOMe) so HCG def is identical (ACG/TCG, GCG+CCG dropped).

## Pipeline audit vs paper (verdict: faithful)
- Trim: `--clip_r1 9 --clip_r2 9 --paired` OK.
- Human depletion: hg38 PE `--non_directional --unmapped` OK (vs paper GRCh37).
- Mouse: SE `--non_directional` on concatenated human-unmapped R1+R2 OK
  (concat-SE ≡ separate-SE; non-directional handles R2; markdup is coord-based
  so identical read names are harmless).
- Dedup before methylation OK.
- Builds: hg38/mm10 vs paper GRCh37/NCBIM37; mm10 matches GEO GRCm38 cov.
- MISSING from paper: high-coverage artifact-region exclusion (not implemented).

## Files created/modified this session
- `acc_list_esc.txt` (NEW) — 32 ESC SRRs.
- `metadata_esc.tsv` (NEW) — SRR→condition→label (2i_1..12, Ser_1..20).
- `make_metadata_esc.sh` (NEW) — builds metadata + subsets QC summary.
- `smallwood_qc_summary_esc.csv` (NEW) — QC summary subset to 32 ESC cells.
- `codes/02.alignment.sh` (MODIFIED) — activated full workflow (hg38 PE depletion
  -> concat unmapped -> mm10 SE), added resume guards. Commit 2164505.
- `codes/05.hcg_nome.sh` (NEW) — coverage2cytosine --nome-seq HCG track, array 1-32.
- `codes/count_hcg_sites.py` (NEW) — per-cell allCpG vs HCG site counts.
- `codes/_test_hcg_one.sh` (NEW) — one-cell sanity test for the HCG step.
- `codes/JOBS.md` (MODIFIED) — added 05.hcg_nome entry.

## Jobs
- 05.hcg_nome (`smwhcg`): script ready, **NOT yet submitted**.
  Submit: `sbatch codes/05.hcg_nome.sh` then run count_hcg_sites.py.

## Next
- Confirm `_test_hcg_one.sh` output (HCG cov rows) then submit 05.hcg_nome.
- Optional: implement high-coverage region exclusion; quantify human-depletion %.
