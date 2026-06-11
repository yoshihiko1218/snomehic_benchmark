#!/usr/bin/env python
"""Count detected CpG (HCG) loci per cell from yap allc-CGN merge files.
Each row in a {cell}.CGN-Merge.allc.tsv.gz = one covered, strand-merged CpG locus.
Usage: count_hcg_loci.py <run_dir> [out_csv]
"""
import sys, gzip, glob, os
import pandas as pd

run_dir = sys.argv[1]
out_csv = sys.argv[2] if len(sys.argv) > 2 else None

files = sorted(glob.glob(os.path.join(run_dir, "**", "*.CGN-Merge.allc.tsv.gz"), recursive=True))
rows = []
for f in files:
    cell = os.path.basename(f).replace(".CGN-Merge.allc.tsv.gz", "")
    n = 0
    mc = 0
    cov = 0
    with gzip.open(f, "rt") as fh:
        for line in fh:
            n += 1
            p = line.rstrip("\n").split("\t")
            # allc cols: chrom pos strand context mc cov methylated
            mc += int(p[4]); cov += int(p[5])
    rows.append({"cell_id": cell, "hcg_loci": n,
                 "mCG": mc, "covCG": cov,
                 "mCGFrac": (mc / cov) if cov else float("nan")})

d = pd.DataFrame(rows)
print(f"run_dir: {run_dir}")
print(f"cells: {len(d)}")
print("HCG (CpG) loci per cell:")
print(f"  median = {d['hcg_loci'].median():,.0f}")
print(f"  mean   = {d['hcg_loci'].mean():,.0f}")
print(f"  min    = {d['hcg_loci'].min():,.0f}")
print(f"  max    = {d['hcg_loci'].max():,.0f}")
print(f"  total  = {d['hcg_loci'].sum():,.0f}")
print(f"global mCG fraction (median per cell) = {d['mCGFrac'].median():.4f}")
if out_csv:
    d.sort_values("hcg_loci", ascending=False).to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")
