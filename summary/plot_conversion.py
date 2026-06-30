#!/usr/bin/env python
"""
Plot per-cell bisulfite-conversion proxy (non-CpG methylation %, = trinuc ACT%)
across methylation methods, as two violin figures matching the qc.ipynb style:
  figures/conversion_chrM_violin.{pdf,png}        (chrM context)
  figures/conversion_chr21chr19_violin.{pdf,png}  (autosome: chr21 human / chr19 mouse)
Reads summary/conversion_percell.csv (built by summary/collect_conversion.py).
"""
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
df = pd.read_csv(f"{BASE}/summary/conversion_percell.csv")

# methylation methods only (no pure-Hi-C nagano/droplethic), low->high coverage order
order = ["scnome", "smallwood", "snmCseq2", "snmCseq3", "snmCAT", "scnomehic"]
palette = {
    "scnome": "#D55E00", "smallwood": "#009E73", "snmCseq2": "#E69F00",
    "snmCseq3": "#56B4E9", "snmCAT": "#7F7F7F", "scnomehic": "#FF0000",
}
labels = {
    "scnome": "scNOMe-seq\n(Pott 2017)", "smallwood": "scWGBS\n(Smallwood 2014)",
    "snmCseq2": "snmC-seq2\n(Luo 2018)", "snmCseq3": "snm3C-seq\n(Liu 2023)",
    "snmCAT": "snmCAT-seq\n(Luo 2022)", "scnomehic": "scNOMe-HiC",
}


def winsorize(d, col):
    out = []
    for g, sub in d.groupby("dataset", sort=False):
        sub = sub.copy()
        q1, q3 = sub[col].quantile(0.25), sub[col].quantile(0.75)
        iqr = q3 - q1
        sub[col] = sub[col].clip(q1 - 1.5 * iqr, q3 + 1.5 * iqr)
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def make(col, title, fname):
    d = df[["dataset", col]].dropna(subset=[col]).copy()
    d = d[d["dataset"].isin(order)]
    counts = d.groupby("dataset").size().to_dict()
    d = winsorize(d, col)
    sns.set_context("paper", font_scale=1.8)
    sns.set_style("whitegrid")
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=d, x="dataset", y=col, order=order, palette=palette,
                   cut=0, inner="box", linewidth=1, ax=ax)
    ax.set_title(title)
    ax.set_ylabel("Non-CpG methylation (%)  [ACT, conversion proxy]")
    ax.set_xlabel("")
    ax.set_xticklabels([f"{labels[o]}\nn={counts.get(o, 0)}" for o in order],
                       rotation=0, fontsize=11)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{BASE}/figures/{fname}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{fname}.{{pdf,png}}  (n per dataset: "
          + ", ".join(f"{o}={counts.get(o,0)}" for o in order) + ")")


make("noncpg_chrM", "Bisulfite non-conversion — chrM", "conversion_chrM_violin")
make("noncpg_auto", "Bisulfite non-conversion — chr21 (human) / chr19 (mouse)",
     "conversion_chr21chr19_violin")
