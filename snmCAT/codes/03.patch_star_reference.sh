#!/bin/bash
# Workaround for a yap 1.6.9 bug in mct mode: the generated Snakefile uses
# {bismark_reference} and {star_reference} in shell rules but yap fails to write
# those values into the Snakefile variable header. Inject them just before the
# CELL_IDS line in every Group Snakefile. Idempotent: inserts if missing, and
# rewrites the value if it points to a different path (e.g. after switching the
# STAR index to the 2.7.10a build that matches yap's STAR version).
set -euo pipefail

BISMARK_REF="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark"
# STAR index MUST be built with the same STAR version yap uses (2.7.10a). The
# 2.7.11b index is unreadable by 2.7.10a (FATAL: unrecognized parameter "genomeType").
STAR_REF="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/star_2.7.10a_gencode.v36_sjdb100"
MAPPING_DIR="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/mapping"

n=0
for sf in "${MAPPING_DIR}"/Group*/Snakefile; do
    # bismark_reference: insert if missing
    if ! grep -qE "^bismark_reference = " "$sf"; then
        sed -i "0,/^CELL_IDS = /s##bismark_reference = '${BISMARK_REF}'\nCELL_IDS = #" "$sf"
    fi
    # star_reference: insert if missing, else overwrite the existing value
    if grep -qE "^star_reference = " "$sf"; then
        sed -i "s#^star_reference = .*#star_reference = '${STAR_REF}'#" "$sf"
    else
        sed -i "0,/^CELL_IDS = /s##star_reference = '${STAR_REF}'\nCELL_IDS = #" "$sf"
    fi
    n=$((n+1))
done
echo "Patched ${n} Snakefiles (bismark_reference + star_reference=${STAR_REF})."
