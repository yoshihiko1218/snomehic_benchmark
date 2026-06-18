#!/usr/bin/env python3
"""
Per-dataset, per-cell summary of HCG and GCH detected loci.

HCG (all 5 methylation datasets): from <dataset>.hcg.txt  (consistent destranded
    + genome GCG removal; Bismark/allcools for the 4, BisSNP NOMe for scnomehic).
GCH (GpC accessibility): NOMe methods ONLY -> scnome (Bismark) + scnomehic (BisSNP);
    from <dataset>.gch.txt. N/A for the bisulfite-only methods (no GpC enzyme).

Writes:
  summary/gch_hcg_counts/<dataset>.summary.txt    per-cell (sample, HCG_n, GCH_n)
  summary/gch_hcg_counts/all_methods.summary.txt  combined (dataset, sample, HCG_n, GCH_n)
  + prints a median table.
"""
import os
import numpy as np
import pandas as pd

OUT = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/summary/gch_hcg_counts"
DATASETS = ["scnome", "smallwood", "snmCseq2", "snmCseq3", "snmCAT", "scnomehic"]
# Datasets to drop LOW outliers from (log-scale lower whisker on HCG, since the
# violins are log-scaled). scnomehic has one failed cell ~50x below the rest.
DROP_LOW_OUTLIERS = {"scnomehic"}
# GCH source per NOMe dataset (GpC, GCG removed, DESTRANDED -- no double counting):
#   scnome    -> Bismark CX cov, GpC destrand + genome GCG removal (scnome.gch.txt)
#   snmCAT    -> YAP allc, GpC destrand + genome GCG removal (snmCAT.gch.txt; mapping_brain)
#   scnomehic -> BisSNP NOMe GCH.6plus2 (destranded by BisSNP) (scnomehic.gch.txt)
GCH_SRC = {"scnome": "scnome.gch.txt", "snmCAT": "snmCAT.gch.txt",
           "scnomehic": "scnomehic.gch.txt"}
# scnome's <ds>.hcg.txt/<ds>.gch.txt are keyed per MATE (K562_01_1, K562_01_2 = R1/R2
# single-end BAMs), which double-counts cells. R1/R2 cover ~the same cytosines, so the
# correct per-cell value is the UNION (~one mate), not the sum. Use the per-cell QC
# summary (already combined: HCG_site_count / GCH_site_count over the 23 finalized cells:
# 12 GM12878 + 11 merged K562). See scnome/SESSION_NOTE_2026-06-09_1.md.
SCNOME_QC = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome/scnome_qc_summary.csv"

combined, med = [], []
for ds in DATASETS:
    if ds == "scnome":
        qc = pd.read_csv(SCNOME_QC)
        d = pd.DataFrame({"sample": qc["CellID"],
                          "HCG_n": qc["HCG_site_count"].astype("Int64"),
                          "GCH_n": qc["GCH_site_count"].astype("Int64")}).sort_values("sample")
        d.to_csv(os.path.join(OUT, f"{ds}.summary.txt"), sep="\t", index=False)
        c = d.copy(); c.insert(0, "dataset", ds)
        combined.append(c[["dataset", "sample", "HCG_n", "GCH_n"]])
        med.append((ds, len(d), int(d["HCG_n"].median()), int(d["GCH_n"].median())))
        continue
    h = pd.read_csv(os.path.join(OUT, f"{ds}.hcg.txt"), sep="\t")[["sample", "HCG_n"]]
    if ds in GCH_SRC:
        g = pd.read_csv(os.path.join(OUT, GCH_SRC[ds]), sep="\t")[["sample", "GCH_n"]]
        d = h.merge(g, on="sample", how="outer")
    else:
        d = h.copy()
        d["GCH_n"] = pd.NA
    d = d.sort_values("sample")
    if ds in DROP_LOW_OUTLIERS:
        lv = np.log10(d["HCG_n"].clip(lower=1))
        q1, q3 = lv.quantile(0.25), lv.quantile(0.75)
        thr = 10 ** (q1 - 1.5 * (q3 - q1))
        n0 = len(d)
        d = d[d["HCG_n"] >= thr]
        if n0 != len(d):
            print(f"  [{ds}] removed {n0 - len(d)} low-outlier cell(s) (HCG < {int(thr):,})")
    d.to_csv(os.path.join(OUT, f"{ds}.summary.txt"), sep="\t", index=False)
    c = d.copy(); c.insert(0, "dataset", ds)
    combined.append(c[["dataset", "sample", "HCG_n", "GCH_n"]])
    gch_med = int(d["GCH_n"].median()) if ds in GCH_SRC else None
    med.append((ds, len(d), int(d["HCG_n"].median()), gch_med))

allc = pd.concat(combined, ignore_index=True)
allc.to_csv(os.path.join(OUT, "all_methods.summary.txt"), sep="\t", index=False)

print(f"{'dataset':10} {'cells':>6} {'median HCG':>12} {'median GCH':>12}")
for ds, n, h, g in med:
    print(f"{ds:10} {n:>6} {h:>12,} {('—' if g is None else format(g, ',')):>12}")
print(f"\nWrote per-dataset <ds>.summary.txt and all_methods.summary.txt ({len(allc)} cells)")
