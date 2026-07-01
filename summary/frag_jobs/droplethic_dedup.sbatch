#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 30:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=dhic_dedup
#SBATCH --output=summary/frag_jobs/logs/%x.%j.out
#SBATCH --error=summary/frag_jobs/logs/%x.%j.err
# droplethic per-barcode dedup: single streaming pass over the position-sorted pairsam.
# A PCR duplicate = same barcode (CB) + same both-end mapping position+strand.
# Emits per-CB: pairs & cis & cis>1kb & trans, BEFORE and AFTER dedup.
source /home/jmj7858/.bashrc 2>/dev/null
conda activate scnomehic 2>/dev/null
set -uo pipefail
B=/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark
cd "$B" || exit 1
PS=droplethic/my_project/03.mapping/SRR27586278_hg38_sorted.pairsam
OUT=summary/frag_counts/droplethic_dedup_percell.tsv

# columns: 1 readID 2 chrom1 3 pos1 4 chrom2 5 pos2 6 strand1 7 strand2 8 pair_type ... 11 CB1
awk -F'\t' '
BEGIN{ print "cell\tpre_pairs\tpre_cis\tpre_cis1kb\tpre_trans\tpost_pairs\tpost_cis\tpost_cis1kb\tpost_trans" > "'"$OUT"'" }
/^#/ { next }
{
  c1=$2; p1=$3; c2=$4; p2=$5; s1=$6; s2=$7; cb=$11;
  if (c1=="!" || c2=="!") next;                     # skip unmapped
  poskey=c1"_"p1"_"c2"_"p2"_"s1"_"s2;
  if (poskey != prevpos) { delete seen; prevpos=poskey }
  isdup = (cb in seen); seen[cb]=1;
  iscis = (c1==c2); d = (p1>p2)?p1-p2:p2-p1;
  pre_p[cb]++; if(iscis){pre_c[cb]++; if(d>1000)pre_k[cb]++} else pre_t[cb]++;
  if(!isdup){ post_p[cb]++; if(iscis){post_c[cb]++; if(d>1000)post_k[cb]++} else post_t[cb]++ }
}
END{
  for(cb in pre_p)
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", cb,
      pre_p[cb], pre_c[cb]+0, pre_k[cb]+0, pre_t[cb]+0,
      post_p[cb]+0, post_c[cb]+0, post_k[cb]+0, post_t[cb]+0 >> "'"$OUT"'"
}' "$PS"
echo "droplethic dedup done -> $OUT ($(wc -l < $OUT) barcodes)"
