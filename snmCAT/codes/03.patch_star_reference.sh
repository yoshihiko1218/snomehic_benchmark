#!/bin/bash
# Workaround for a yap 1.6.9 bug in mct mode: the generated Snakefile uses
# {bismark_reference} and {star_reference} in shell rules but yap fails to write
# those values into the Snakefile variable header. Inject them just before the
# CELL_IDS line in every Group Snakefile. Idempotent.
set -euo pipefail

BISMARK_REF="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark"
STAR_REF="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/star_2.7.11b_gencode.v36_sjdb150"
MAPPING_DIR="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/mapping"

n=0
for sf in "${MAPPING_DIR}"/Group*/Snakefile; do
    patched=0
    if ! grep -qE "^bismark_reference = " "$sf"; then
        sed -i "0,/^CELL_IDS = /s##bismark_reference = '${BISMARK_REF}'\nCELL_IDS = #" "$sf"
        patched=1
    fi
    if ! grep -qE "^star_reference = " "$sf"; then
        sed -i "0,/^CELL_IDS = /s##star_reference = '${STAR_REF}'\nCELL_IDS = #" "$sf"
        patched=1
    fi
    n=$((n+patched))
done
echo "Patched ${n} Snakefiles (bismark_reference + star_reference)."
