#!/bin/bash
set -eo pipefail
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark

git add snmCseq3/codes/comparison/PIPELINE_COMPARISON.md
git commit -m "Fix --nofw/--norc description: flags control query orientation, not reference" \
  -m "" \
  -m "Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>"
GIT_ASKPASS="" git push origin main
echo "Done"
