#!/usr/bin/env python3
"""
scnome per-cell UNIQUE HCG and GCH detected-loci, MERGED across the two R1/R2 mates.

scnome was aligned as two single-end mates (_1/_2), so its covs are per-mate and the
canonical count_hcg_destranded / count_gch count each mate separately (double rows).
The correct per-cell value is the UNION of the two mates' destranded loci sets.

Same input + detection as the canonical scripts (one CX cov per mate):
  HCG: CpG destranded, GCG removed     (== count_hcg_destranded.count_destranded)
  GCH: GpC destranded, GCG removed     (== count_gch.gch_destranded)
23 cells: GM SRR3729642-3729653 + merged K562_01-11 (raw K562 SRR>653 + controls dropped).

Output: summary/gch_hcg_counts/scnome.loci_percell.txt  (sample, HCG_n, GCH_n)
"""
import os, glob, gzip, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/summary")
import hcg_lib
import pandas as pd

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
METHY = os.path.join(BASE, "scnome/05.methy")
OUT = os.path.join(BASE, "summary/gch_hcg_counts/scnome.loci_percell.txt")
SUF = ".rmdup.bismark.cov.gz"


def cell_loci(args):
    cell, paths = args
    seqs = hcg_lib._GENOMES["hg38"]
    hcg, gch = set(), set()
    for path in paths:                         # union across the 2 mates
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.split("\t", 2)
                chrom = f[0]
                seq = seqs.get(chrom)
                if seq is None:
                    continue
                i = int(f[1]) - 1
                if i < 2 or i + 2 >= len(seq):
                    continue
                b = seq[i]
                # --- HCG: CpG dinucleotide, destranded, GCG removed ---
                if b == "C" and seq[i + 1] == "G":
                    cpg = i
                elif b == "G" and seq[i - 1] == "C":
                    cpg = i - 1
                else:
                    cpg = -1
                if cpg >= 1 and seq[cpg - 1] != "G":
                    hcg.add((chrom, cpg))
                # --- GCH: GpC dinucleotide, destranded, GCG removed ---
                if b == "C" and seq[i - 1] == "G":
                    g = i - 1; is_gcg = seq[i + 1] == "G"
                elif b == "G" and seq[i + 1] == "C":
                    g = i; is_gcg = seq[i + 2] == "G"
                else:
                    g = -1; is_gcg = True
                if g >= 0 and not is_gcg:
                    gch.add((chrom, g))
    return cell, len(hcg), len(gch)


def build_cells():
    bycell = defaultdict(list)
    for f in glob.glob(os.path.join(METHY, "*_[12]" + SUF)):
        s = os.path.basename(f)[:-len(SUF)]
        cell = s.rsplit("_", 1)[0]
        if cell.startswith("SRR"):
            try:
                if int(cell[3:]) > 3729653:
                    continue
            except ValueError:
                pass
        bycell[cell].append(f)
    return {c: sorted(p) for c, p in bycell.items()}


def main():
    hcg_lib.load_genomes({"hg38": hcg_lib.HG38})
    bycell = build_cells()
    print(f"[scnome per-cell HCG+GCH] {len(bycell)} cells", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for cell, h, g in ex.map(cell_loci, list(bycell.items())):
            rows.append((cell, h, g))
            print(f"  {cell}: HCG={h:,} GCH={g:,}", flush=True)
    df = pd.DataFrame(rows, columns=["sample", "HCG_n", "GCH_n"]).sort_values("sample")
    df.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} ({len(df)} cells); median HCG={df['HCG_n'].median():.0f} "
          f"GCH={df['GCH_n'].median():.0f}", flush=True)


if __name__ == "__main__":
    main()
