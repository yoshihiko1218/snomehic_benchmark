#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=keepalive_touch
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scratch_keepalive/logs/keepalive_%j.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scratch_keepalive/logs/keepalive_%j.out

# ---------------------------------------------------------------------------
# Quest b1042 scratch keep-alive.
# Refreshes atime+mtime on the whole yoshii/ tree so the 30-day inactivity
# purge never flags it, then resubmits itself to run again in 7 days.
# Stop the chain by: scancel'ing the pending keepalive_touch job (squeue -u $USER).
# ---------------------------------------------------------------------------

set -uo pipefail

SCRIPT="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scratch_keepalive/keepalive_touch.sh"
TARGET="/gpfs/projects/b1042/epifluidlab/yoshii/"
INTERVAL_DAYS=7

echo "=== keepalive run start: $(date) on $(hostname) (job ${SLURM_JOB_ID:-NA}) ==="

# 1) Refresh access + modify times across the whole tree.
find "$TARGET" -exec touch -a -m {} + 2> >(grep -v '^$' >&2)
rc=$?
nfiles=$(find "$TARGET" | wc -l)
echo "touched tree (exit $rc); current file/dir count: $nfiles"

# 2) Resubmit self to run again in INTERVAL_DAYS, so this is perpetual.
NEXT=$(date -d "+${INTERVAL_DAYS} days" +%Y-%m-%dT%H:%M:%S)
newjob=$(sbatch --begin="$NEXT" "$SCRIPT" 2>&1)
echo "resubmit for $NEXT -> $newjob"

echo "=== keepalive run end: $(date) ==="
