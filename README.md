# scnomehic_paper / benchmark

Cross-technology benchmark of single-cell methylation / accessibility / Hi-C
methods. Each top-level folder is one dataset (raw → aligned → QC), and `summary/`
holds the cross-method QC notebook and comparison outputs.

## Layout

### Dataset folders (pipeline working dirs — internal layout left as-is)
| Folder | Technology / source | Size |
|--------|---------------------|------|
| `scnomehic/`     | scNOMe-HiC (this method); GM12878 + others | 707G |
| `nagano/`        | Nagano scHi-C | 102G |
| `droplethic/`    | Droplet Hi-C | 1.6T |
| `scnome/`        | scNOMe | 839G |
| `smallwood/`     | Smallwood scBS-seq (ESC-only subset) | 880G |
| `snmCAT/`        | snmC2T-seq / NOMe (brain) | 66G |
| `snmCseq2/`      | snmC-seq2 | 1.1T |
| `snmCseq3/`      | snm3C-seq | 549G |
| `methylhic/`     | Methyl-HiC | 785G |
| `methylhic_new/` | Methyl-HiC (newer batch) | 65G |

Each dataset folder has a `codes/` (pipeline scripts), alignment dir(s), `logs/`,
and most have `qc_stats/`. See each folder's `FILES.md` where present. Internal
folder names are intentionally NOT standardized — Snakefiles and `.snakemake/`
resume state hardcode these paths (see `REORG_PLAN.md`).

### Support folders
- `summary/` — `qc.ipynb` (master cross-method QC notebook) + comparison scripts/outputs.
  Real inputs the notebook reads live in `summary/trinuc/`, `summary/hic_cache/`,
  `summary/gch_hcg_counts/`. Loose helper scripts stay here (coupled via `import hcg_lib`
  and `submit_*.sh`); one-off notebook-patch scripts are archived in `summary/_archive/`.
- `figures/` — published QC figures (PDF/PNG) written by `qc.ipynb`.
- `tools/` — repo-level helper scripts (git add/commit/push wrappers). Run from the
  benchmark root, e.g. `bash tools/commit_all.sh`.
- `notes/` — dated `SESSION_NOTE_*.md` work logs.
- `scratch_keepalive/` — weekly Slurm job that touches the b1042 scratch tree to
  avoid the 30-day inactivity purge (see its README).
- `logs/` — misc job logs.

### Root-level shared tools (kept at root; referenced by absolute path from pipelines)
- `bam_summary_universal.py` — BAM → alignment/Hi-C summary (called by `scnome/codes/*.sh`).
- `build_trinuc_summary.py` + `parse_bissnp_trinuc.py` — trinucleotide methylation summary.
- `trinuc_qc_summary.csv` — output of the above.
- `PROJECT_CONTEXT.md` — project background.

## Conventions
- New session notes go in `notes/`.
- Per-folder `FILES.md` documents what each file is and how it was generated.
- Submitted jobs are tracked in each folder's `JOBS.md`.
