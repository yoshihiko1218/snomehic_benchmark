#!/usr/bin/env python
"""
Collect per-cell bisulfite-conversion proxy (non-CpG methylation %, = trinuc ACT%)
for the 6 methylation datasets, for two genomic contexts:
  - chrM   (mitochondrial, naked DNA -> residual mC = non-conversion proxy)
  - autosome (chr21 for human / chr19 for mouse) -> autosomal non-CpG proxy
and write a tidy per-cell table summary/conversion_percell.csv with columns:
  dataset, cell, noncpg_chrM, noncpg_auto, auto_chrom

Per-cell, NEVER per-mate. Cell sets match the other QC panels:
  scnome 23 | smallwood 51 | snmCseq2 96 (mm10) | snmCseq3 98 | scnomehic 187 | snmCAT 99

Sources (per-cell unless noted):
  scnome    : scnome_qc_summary.csv  (chrM_noncpg, chr21_noncpg)            [Bismark]
  smallwood : smallwood_qc_summary.csv (chrM_noncpg, chr19_noncpg)          [Bismark]
  snmCseq2  : chrM from snmcseq2_qc_summary.csv (mm10 cells); chr19 from the
              PER-MATE trinuc/snmCseq2.chr19.txt averaged per cell           [Bismark]
  snmCseq3  : trinuc/snmCseq3.{chrM,chr21}.txt                              [bhmem]
  scnomehic : external gm.{chrM,chr21}.txt.gz, restricted to gm_passed.txt   [BisSNP]
  snmCAT    : trinuc/snmCAT.{chrM,chr21}.txt (computed from YAP allc),
              restricted to the 99 cells in gch_hcg_counts/all_methods.summary.txt [YAP allc]
NOTE: snmCAT is the one position-level (allc Sigma mc / Sigma cov) source; all others are
read-level ACT% from BAM trinuc. Same ACT/ACG/GCT proxy, slightly different weighting.
"""
import os
import numpy as np
import pandas as pd

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
TRI = os.path.join(BASE, "summary/trinuc")
EXT = "/home/jmj7858/epifluidlab/scnomehic_paper/s1/QC"
rows = []


def add(ds, cell, chrM, auto, auto_chrom):
    rows.append({"dataset": ds, "cell": str(cell),
                 "noncpg_chrM": chrM, "noncpg_auto": auto, "auto_chrom": auto_chrom})


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


# --- scnome (23) ---
q = pd.read_csv(f"{BASE}/scnome/scnome_qc_summary.csv")
for _, r in q.iterrows():
    add("scnome", r["CellID"], num(r.get("chrM_noncpg")), num(r.get("chr21_noncpg")), "chr21")

# --- smallwood (51) ---
q = pd.read_csv(f"{BASE}/smallwood/smallwood_qc_summary.csv")
for _, r in q.iterrows():
    add("smallwood", r["CellID"], num(r.get("chrM_noncpg")), num(r.get("chr19_noncpg")), "chr19")

# --- snmCseq2 (96 mm10): chrM from qc_summary, chr19 from per-mate trinuc averaged ---
mm10 = set(pd.read_csv(f"{BASE}/snmCseq2/codes/cells_mm10.txt", header=None)[0].astype(str))
q = pd.read_csv(f"{BASE}/snmCseq2/snmcseq2_qc_summary.csv")
q = q[q["CellID"].astype(str).isin(mm10)]
chrm_map = {str(r["CellID"]): num(r.get("chrM_noncpg")) for _, r in q.iterrows()}
t = pd.read_csv(f"{TRI}/snmCseq2.chr19.txt", sep="\t")
t["cell"] = t["sample"].astype(str).str.replace(r"_[12]$", "", regex=True)
auto_map = t.groupby("cell")["noncpg"].mean().to_dict()   # average R1/R2 mates
for cell in sorted(mm10):
    add("snmCseq2", cell, chrm_map.get(cell, np.nan), auto_map.get(cell, np.nan), "chr19")

# --- snmCseq3 (98): both from trinuc ---
cm = pd.read_csv(f"{TRI}/snmCseq3.chrM.txt", sep="\t").set_index("sample")["noncpg"].to_dict()
au = pd.read_csv(f"{TRI}/snmCseq3.chr21.txt", sep="\t").set_index("sample")["noncpg"].to_dict()
for cell in au:   # chr21 file defines the 98-cell set
    add("snmCseq3", cell, num(cm.get(cell)), num(au.get(cell)), "chr19")  # mm10 -> chr19 biology

# --- scnomehic (187 passing) ---
passed = set(pd.read_csv(f"{EXT}/gm_passed.txt", header=None)[0].astype(str))
cm = pd.read_csv(f"{EXT}/chrM_gch_hcg/gm.chrM.txt.gz", sep="\t").set_index("sample")["noncpg"].to_dict()
au = pd.read_csv(f"{EXT}/chrM_gch_hcg/gm.chr21.txt.gz", sep="\t").set_index("sample")["noncpg"].to_dict()
for cell in passed:
    add("scnomehic", cell, num(cm.get(cell)), num(au.get(cell)), "chr21")

# --- snmCAT (99): restrict to the loci-panel cell set ---
am = pd.read_csv(f"{BASE}/summary/gch_hcg_counts/all_methods.summary.txt", sep="\t")
snmcat_cells = set(am[am["dataset"].str.lower() == "snmcat"]["sample"].astype(str))
cm = pd.read_csv(f"{TRI}/snmCAT.chrM.txt", sep="\t").set_index("sample")["noncpg"].to_dict()
au = pd.read_csv(f"{TRI}/snmCAT.chr21.txt", sep="\t").set_index("sample")["noncpg"].to_dict()
use = snmcat_cells if snmcat_cells else set(au.keys())
for cell in sorted(use):
    add("snmCAT", cell, num(cm.get(cell)), num(au.get(cell)), "chr21")

df = pd.DataFrame(rows)
out = f"{BASE}/summary/conversion_percell.csv"
df.to_csv(out, index=False)

print(f"wrote {out}  ({len(df)} cell-rows)")
print("\nper-dataset cell counts & median non-CpG %:")
g = df.groupby("dataset").agg(
    n=("cell", "size"),
    n_chrM=("noncpg_chrM", lambda s: s.notna().sum()),
    med_chrM=("noncpg_chrM", "median"),
    n_auto=("noncpg_auto", lambda s: s.notna().sum()),
    med_auto=("noncpg_auto", "median"),
)
print(g.round(3).to_string())
