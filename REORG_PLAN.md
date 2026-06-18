# Benchmark reorganization plan

Goal: make `benchmark/` more organized. Investigation found heavy path coupling
(qc.ipynb relative paths; hundreds of auto-generated Snakefiles with absolute
paths + live .snakemake state; codes/*.sh BASE vars). This plan is tiered by risk.

## Current pain points
- **Root clutter:** git helpers (`add_all.sh`, `commit_all.sh`, `push_*.sh`, `add_specific.sh`),
  utility py (`bam_summary_universal.py`, `build_trinuc_summary.py`, `parse_bissnp_trinuc.py`),
  loose `trinuc_qc_summary.csv`, 4 old `SESSION_NOTE_*.md`.
- **summary/ clutter:** 24 helper `.py` (incl. one-off `_fix_*.py`/`_add_*.py`/`*_nb.py`
  notebook-patchers), 13 `.sh`, 11 output CSVs, 9 PDFs, `.bak` notebooks. Real subfolders
  (`trinuc/`, `hic_cache/`, `gch_hcg_counts/`) are referenced by qc.ipynb and must stay.
- **Inconsistent dataset layouts:** numbered (scnome, smallwood, snmCseq2, methylhic) vs flat
  (nagano, snmCAT, scnomehic, methylhic_new) vs chaotic (snmCseq3: 7+ `04.bhmem_bam_*` variants).
- **Docs gaps:** nagano, droplethic, methylhic, methylhic_new have no README/FILES.md.

## Tier 1 — tidy + document (LOW risk, recommended, no pipeline breakage)
Dataset folders and ALL their internals stay exactly where they are (so Snakefiles,
.snakemake state, and qc.ipynb relative paths are untouched).
- Root: `tools/` for the 6 git helper `.sh` + 3 utility `.py`; move `trinuc_qc_summary.csv`
  into `summary/`; `notes/` for the 4 root `SESSION_NOTE_*.md`; add top-level `README.md` index.
- summary/: `summary/scripts/` for helper `.py`/`.sh`; `summary/scripts/_archive/` for one-off
  `_fix_*.py`/`_add_*.py`/`*_nb.py`; `summary/_backups/` for `.bak` notebooks. Keep `qc.ipynb`,
  `trinuc/`, `hic_cache/`, `gch_hcg_counts/` in place. Leave the big output CSVs/PDFs in summary/.
- Add a short `FILES.md` to each undocumented dataset folder.
- Rewrite ONLY the handful of internal path refs broken by the above moves (verify with grep).

## Tier 2 — standardize dataset internals (MEDIUM-HIGH risk, optional)
Rename each dataset to a uniform stage scheme (e.g. 01.fastq/02.trim/04.align/06.methy/qc_stats),
consolidate snmCseq3 alignment variants. Requires sed-rewriting hundreds of Snakefiles + codes
scripts + qc.ipynb refs, and invalidates .snakemake resume state. Must be done per-dataset with
testing. High effort, breakage-prone. Recommendation: document structure instead, or defer.

## Tier 3 — full restructure (HIGHEST risk, not recommended)
Move all datasets under `datasets/` and rewrite every relative + absolute ref. Aesthetic gain
only; large blast radius. Not recommended.

## Decision
User chose **Tier 1 only** (2026-06-18). Tier 2 & 3 skipped (pipeline-path risk).

## Tier 1 — DONE (2026-06-18)
- `tools/` ← 6 git helper scripts (run from repo root, e.g. `bash tools/commit_all.sh`).
- `notes/` ← 4 old root `SESSION_NOTE_*.md`.
- `summary/_archive/` ← 8 one-off notebook-patcher scripts (verified unreferenced).
- `summary/_backups/` ← `.bak` notebooks (gitignored).
- Added root `README.md`; `FILES.md` to nagano/droplethic/methylhic/methylhic_new.
- Left coupled tools in place (`bam_summary_universal.py`, `build_trinuc_summary.py`,
  `parse_bissnp_trinuc.py`, `trinuc_qc_summary.csv`, all working `summary/` scripts).
- No dataset-folder internals touched; no Snakefiles/.snakemake modified.
