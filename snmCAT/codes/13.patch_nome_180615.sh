#!/bin/bash
# Patch mapping_180615 Snakefiles: inject bismark_reference + star_reference (yap
# omits them) and set nome_flag_str='--nome' (yap leaves it blank). Idempotent.
set -euo pipefail
BISMARK_REF="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark"
STAR_REF="/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/star_2.7.10a_gencode.v36_sjdb100"
MAPPING_DIR="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/mapping_180615"
n=0
for sf in "${MAPPING_DIR}"/Group*/Snakefile; do
    grep -qE "^bismark_reference = " "$sf" || sed -i "0,/^CELL_IDS = /s##bismark_reference = '${BISMARK_REF}'\nCELL_IDS = #" "$sf"
    grep -qE "^star_reference = " "$sf"    || sed -i "0,/^CELL_IDS = /s##star_reference = '${STAR_REF}'\nCELL_IDS = #" "$sf"
    sed -i "s#^nome_flag_str = .*#nome_flag_str = '--nome'#" "$sf"
    n=$((n+1))
done
echo "Patched ${n} Snakefiles (bismark+star refs, nome_flag_str='--nome')."
