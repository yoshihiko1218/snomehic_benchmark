#!/usr/bin/env python3
"""
scnome destranded HCG, MERGED per cell across the two R1/R2 mates.

count_hcg_destranded.py counts each mate (_1/_2) cov separately, so a cell's HCG
appears twice; you can't just sum them (R1/R2 cover ~the same CpGs). The correct
per-cell value is the UNION of the two mates' destranded HCG CpG sets.

This reads BOTH mate covs for each of the 23 finalized cells (GM SRR3729642-3729653
+ merged K562_01-11; raw K562 SRRs >3729653 and controls dropped), accumulates the
HCG/GCG CpG sets across both mates, then counts unique CpGs.

Destranding + GCG removal logic is identical to count_hcg_destranded.count_destranded
(Bismark CpG cov; collapse +/- to the CpG + strand C; GCG if ref[cpg-1]==G else HCG).

Output: summary/gch_hcg_counts/scnome.hcg_percell_destranded.txt  (sample, HCG_n, GCG_n)
"""
import os, glob, gzip, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/summary")
import hcg_lib
import pandas as pd

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
METHY = os.path.join(BASE, "scnome/05.methy")
OUT = os.path.join(BASE, "summary/gch_hcg_counts/scnome.hcg_percell_destranded.txt")
SUF = ".rmdup.bismark.cov.gz"


def cell_hcg(args):
    cell, paths = args
    seqs = hcg_lib._GENOMES["hg38"]
    hcg, gcg = set(), set()
    for path in paths:                      # union across the 2 mates
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.split("\t", 2)
                seq = seqs.get(f[0])
                if seq is None:
                    continue
                i = int(f[1]) - 1
                if i < 2 or i + 1 >= len(seq):
                    continue
                b = seq[i]
                if b == "C" and seq[i + 1] == "G":
                    cpg = i
                elif b == "G" and seq[i - 1] == "C":
                    cpg = i - 1
                else:
                    continue
                if cpg < 1:
                    continue
                (gcg if seq[cpg - 1] == "G" else hcg).add((f[0], cpg))
    return cell, len(hcg), len(gcg)


def build_cells():
    bycell = defaultdict(list)
    for f in glob.glob(os.path.join(METHY, "*_[12]" + SUF)):
        s = os.path.basename(f)[:-len(SUF)]     # e.g. K562_01_1 / SRR3729642_1
        cell = s.rsplit("_", 1)[0]              # K562_01 / SRR3729642
        if cell.startswith("SRR"):
            try:
                if int(cell[3:]) > 3729653:     # raw K562 SRRs + controls -> skip
                    continue
            except ValueError:
                pass
        bycell[cell].append(f)
    return {c: sorted(p) for c, p in bycell.items()}


def main():
    hcg_lib.load_genomes({"hg38": hcg_lib.HG38})
    bycell = build_cells()
    print(f"[scnome per-cell merge] {len(bycell)} cells "
          f"({sum(len(v) for v in bycell.values())} mate covs)", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for cell, h, g in ex.map(cell_hcg, list(bycell.items())):
            rows.append((cell, h, g))
            print(f"  {cell}: HCG={h:,} GCG={g:,}", flush=True)
    df = pd.DataFrame(rows, columns=["sample", "HCG_n", "GCG_n"]).sort_values("sample")
    df.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} ({len(df)} cells); median HCG={df['HCG_n'].median():.0f}", flush=True)


if __name__ == "__main__":
    main()
