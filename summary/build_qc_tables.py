#!/usr/bin/env python
"""
Assemble the CONSISTENT per-cell QC tables that qc.ipynb plots, from the harmonized
sources. Writes:
  summary/qc_tables/align.csv  -> dataset, cell, uniq_preDedup, uniq_postDedup,
                                   rate_preDedup, rate_postDedup   (metrics 2 & 3)
  summary/qc_tables/hic.csv    -> dataset, cell, cis_n, cis_gt1kb, trans_ratio,
                                   cis_per_million                 (metrics 4,5,6; MapQ30)
Conversion (metric 1) and loci already live in summary/conversion_percell.csv and
summary/gch_hcg_counts/all_methods.summary.txt.
"""
import glob, gzip, os
import numpy as np
import pandas as pd

B = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
GM = "/home/jmj7858/epifluidlab/sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake"
OUT = f"{B}/summary/qc_tables"
os.makedirs(OUT, exist_ok=True)

# ---------- metrics 2 & 3 ----------
fc = pd.read_csv(f"{B}/summary/frag_counts_all.tsv", sep="\t")
fc[["dataset", "cell", "uniq_preDedup", "uniq_postDedup",
    "rate_preDedup", "rate_postDedup"]].to_csv(f"{OUT}/align.csv", index=False)
print(f"align.csv: {len(fc)} cells")

# ---------- metrics 4,5,6 (Hi-C contacts, MapQ30-consistent) ----------
def parse_summary(path, opener=open):
    d = {}
    with opener(path, "rt") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1] not in ("", None):
                try:
                    d[p[0].rstrip(":")] = float(p[1])
                except ValueError:
                    pass
    return d

rows = []
def add(ds, cell, cis, cis1kb, trans, mapped_mapq30):
    denom = cis + trans
    tr = trans / denom if denom and denom > 0 else np.nan
    cpm = cis1kb / (mapped_mapq30 / 1e6) if mapped_mapq30 and mapped_mapq30 > 0 else np.nan
    rows.append(dict(dataset=ds, cell=str(cell), cis_n=cis, cis_gt1kb=cis1kb,
                     trans_ratio=tr, cis_per_million=cpm))

# nagano
for f in glob.glob(f"{B}/nagano/alignment/*.summary.txt"):
    d = parse_summary(f)
    add("nagano", os.path.basename(f)[:-len(".summary.txt")],
        d.get("UniqMappedMapQ30NoPcrCis", np.nan), d.get("UniqMappedMapQ30NoPcrCisMore1kb", np.nan),
        d.get("UniqMappedMapQ30NoPcrTrans", np.nan), d.get("UniqMappedMapQ30", np.nan))
# scnomehic (187 passed)
passed = set(pd.read_csv("/home/jmj7858/epifluidlab/scnomehic_paper/s1/QC/gm_passed.txt", header=None)[0].astype(str))
for c in passed:
    f = f"{GM}/{c}.summary.txt.gz"
    if os.path.exists(f):
        d = parse_summary(f, gzip.open)
        add("scnomehic", c, d.get("UniqMappedMapQ30NoPcrCis", np.nan), d.get("UniqMappedMapQ30NoPcrCisMore1kb", np.nan),
            d.get("UniqMappedMapQ30NoPcrTrans", np.nan), d.get("UniqMappedMapQ30", np.nan))
# droplethic (valid.tsv) -- NOTE: pre-dedup (dedup never run); flagged in analysis_pipeline.md
dt = pd.read_csv(f"{B}/droplethic/my_project/SRR27586278_hg38.per_cell_qc.valid.tsv", sep="\t")
cc = {c.lower(): c for c in dt.columns}
gc = lambda *n: next((cc[x.lower()] for x in n if x.lower() in cc), None)
c_cis, c_1kb, c_tr, c_mq = gc("UniqMappedMapQ30NoPcrCis"), gc("UniqMappedMapQ30NoPcrCisMore1kb", "UniqMappedMapQ30NoDup_Cis1kb"), gc("UniqMappedMapQ30NoPcrTrans"), gc("UniqMappedMapQ30")
cbcol = gc("CB", "cell", "barcode")
for _, r in dt.iterrows():
    add("droplethic", r[cbcol], r[c_cis], r[c_1kb], r[c_tr], r.get(c_mq, np.nan))
# snmCseq3 MapQ30 (regenerated counts.txt)
for f in glob.glob(f"{B}/summary/frag_counts/snmCseq3_contacts_q30/*.q30.counts.txt"):
    d = dict(l.strip().split(",") for l in open(f) if "," in l)
    short, lng, tr = float(d["CisShortContact"]), float(d["CisLongContact"]), float(d["TransContact"])
    cell = os.path.basename(f)[:-len(".q30.counts.txt")]
    add("snmCseq3", cell, short + lng, lng, tr, short + lng + tr)  # normalize by total contacts

hic = pd.DataFrame(rows)
hic.to_csv(f"{OUT}/hic.csv", index=False)
print(f"hic.csv: {len(hic)} cells")
print(hic.groupby("dataset").agg(n=("cell", "size"), med_cis1kb=("cis_gt1kb", "median"),
      med_trans=("trans_ratio", "median")).round(4).to_string())
