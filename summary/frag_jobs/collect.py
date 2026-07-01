#!/usr/bin/env python
"""Collect per-cell fragment-QC tsvs into one table: summary/frag_counts_all.tsv
Both before-dedup (preDedup) and after-dedup (postDedup) versions.
Columns: dataset, cell, uniq_preDedup, mapq30_preDedup, rate_preDedup,
                        uniq_postDedup, mapq30_postDedup, rate_postDedup"""
import glob, os, pandas as pd

B = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
cols = ["uniq_preDedup", "mapq30_preDedup", "rate_preDedup",
        "uniq_postDedup", "mapq30_postDedup", "rate_postDedup"]
rows = []
for ds in ["nagano", "smallwood", "scnome", "snmCseq3", "scnomehic", "snmCseq2", "snmCAT"]:
    for f in glob.glob(f"{B}/summary/frag_counts/{ds}/*.tsv"):
        parts = open(f).read().strip().split("\t")
        if len(parts) != 7:
            print("skip malformed", os.path.basename(f)); continue
        c = parts[0]
        d = {"dataset": ds, "cell": c}
        d.update({k: (int(v) if i < 2 or i in (3, 4) else float(v))
                  for i, (k, v) in enumerate(zip(cols, parts[1:]))})
        rows.append(d)
dp = f"{B}/summary/frag_counts/droplethic_percell.tsv"
if os.path.exists(dp):
    d = pd.read_csv(dp, sep="\t")
    for _, x in d.iterrows():
        rows.append({"dataset": "droplethic", "cell": str(x["cell"]),
                     **{k: x[k] for k in cols}})

df = pd.DataFrame(rows)
out = f"{B}/summary/frag_counts_all.tsv"
df.to_csv(out, sep="\t", index=False)
print(f"wrote {out}: {len(df)} cells")
print(df.groupby("dataset").agg(
    n=("cell", "size"),
    cnt_pre=("uniq_preDedup", "median"), cnt_post=("uniq_postDedup", "median"),
    rate_pre=("rate_preDedup", "median"), rate_post=("rate_postDedup", "median"),
).round(1).to_string())
