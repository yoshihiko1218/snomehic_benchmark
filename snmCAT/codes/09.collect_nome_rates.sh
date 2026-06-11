#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=nome_rates
#SBATCH --output=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/09.nome_rates/collect.out
#SBATCH --error=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT/codes/logs/09.nome_rates/collect.err

# Aggregate per-cell HCG/GCH/HCH site counts + methylation rates from the
# --nome ALLC files (yap context = upstream(1) + C + downstr(2) = 4 chars).
#   HCG = c[1]in{A,C,T} & c[3]=='G'   (H C G ...)   -> CpG methylation
#   GCH = c[1]=='G'     & c[3]in{A,C,T} (G C H ...) -> NOMe accessibility (if any)
#   HCH = c[1]in{A,C,T} & c[3]in{A,C,T}             -> background
#   (1-indexed awk substr: pos1=upstream, pos2='C', pos3=downstr1)
set -o pipefail   # no -e/-u/bashrc: keep going past any single unreadable file (logged)
BASE="/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT"
OUT="${BASE}/mapping_nome/stats/hcg_gch_nome.tsv"
mkdir -p "${BASE}/mapping_nome/stats" "${BASE}/codes/logs/09.nome_rates"

printf "cell_id\tHCG_site_count\tGCH_site_count\tHCH_site_count\tHCG_mc_rate\tGCH_mc_rate\tHCH_mc_rate\n" > "${OUT}"

for allc in "${BASE}"/mapping_nome/Group*/allc/*.allc.tsv.gz; do
    cell=$(basename "${allc}" .allc.tsv.gz)
    if ! gzip -t "${allc}" 2>/dev/null; then
        echo "WARN: unreadable, skipping ${allc}" >&2
        continue
    fi
    zcat "${allc}" | awk -F'\t' -v cell="${cell}" '
        { c=$4; mc=$5; cov=$6; u=substr(c,1,1); d=substr(c,3,1);
          if ((u=="A"||u=="C"||u=="T") && d=="G")      { hn++; hmc+=mc; hcov+=cov }
          else if (u=="G" && (d=="A"||d=="C"||d=="T")) { gn++; gmc+=mc; gcov+=cov }
          else if ((u=="A"||u=="C"||u=="T") && (d=="A"||d=="C"||d=="T")) { bn++; bmc+=mc; bcov+=cov } }
        END { printf "%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n", cell, hn, gn, bn,
              (hcov?100*hmc/hcov:0), (gcov?100*gmc/gcov:0), (bcov?100*bmc/bcov:0) }' >> "${OUT}"
done

n=$(($(wc -l < "${OUT}") - 1))
echo "Wrote ${OUT} with ${n} cells"
echo "=== medians (nome) ==="
tail -n +2 "${OUT}" | awk -F'\t' '{h[NR]=$5;g[NR]=$6;b[NR]=$7}
  END{asort(h);asort(g);asort(b);m=int(NR/2);
      printf "HCG_mc_rate ~ %.1f%% ; GCH_mc_rate ~ %.2f%% ; HCH(bg) ~ %.2f%%\n", h[m], g[m], b[m]}'
