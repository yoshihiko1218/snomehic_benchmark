#!/bin/bash
# Manual cleanup for the 2026-06-18 QC/loci-fix session.
# Review and run the sections you want:  bash cleanup_session.sh
# (Nothing here is required for the notebook to work; it's just removing cruft.)
set -u
B=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
cd "$B" || exit 1

echo "=== 1) temp notebook-execution copy (safe; recreated by run_qc_notebook.sh) ==="
ls -la _tmp_qc_run.ipynb 2>/dev/null && rm -f _tmp_qc_run.ipynb && echo "  removed _tmp_qc_run.ipynb"

echo "=== 2) my edit backups (gitignored; only needed if you want to revert an edit) ==="
ls -la summary/_backups/ 2>/dev/null
# rm -f summary/_backups/qc.ipynb.bak_*          # uncomment to delete ALL edit backups
# rm -f summary/gch_hcg_counts/all_methods.summary.txt.bak_before_scnome_fix

echo "=== 3) stale hic_df cache (rebuilds on next full notebook run) ==="
ls -la summary/hic_cache/hic_df.csv 2>/dev/null
# rm -f summary/hic_cache/hic_df.csv             # uncomment to force a fresh HiC rebuild

echo "=== 4) superseded HCG-ONLY intermediate (replaced by scnome_loci_percell.*) ==="
echo "    These are committed, so use git rm (then commit) if you want them gone:"
echo "      git rm summary/scnome_hcg_percell.py summary/submit_scnome_hcg_percell.sh \\"
echo "             summary/gch_hcg_counts/scnome.hcg_percell_destranded.txt"
echo "      git commit -m 'remove superseded HCG-only scnome merge intermediates'"

echo "=== done. Nothing was deleted except _tmp_qc_run.ipynb; uncomment the rm lines you want. ==="
