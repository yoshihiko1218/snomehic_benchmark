#!/usr/bin/env python
"""Cross-technology distribution sanity-check of the (now-consistent) metrics.
Prints n / median / Q1 / Q3 / min / max per technology for each metric so we can
eyeball whether the distributions are sensible and comparable."""
import glob, gzip, os
import numpy as np
import pandas as pd

B = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
GM = "/home/jmj7858/epifluidlab/sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake"


def q(s):
    s = pd.to_numeric(pd.Series(list(s)), errors="coerce").dropna()
    if len(s) == 0:
        return None
    return dict(n=len(s), med=s.median(), q1=s.quantile(.25), q3=s.quantile(.75),
                lo=s.min(), hi=s.max())


def show(title, per_tech):
    print(f"\n### {title}")
    print(f"{'tech':12s} {'n':>5s} {'median':>12s} {'Q1':>12s} {'Q3':>12s} {'min':>12s} {'max':>12s}")
    for t, s in per_tech.items():
        d = q(s)
        if d is None:
            print(f"{t:12s}  (no data)"); continue
        f = lambda x: f"{x:12.3f}" if abs(x) < 1000 else f"{x:12.0f}"
        print(f"{t:12s} {d['n']:5d} {f(d['med'])} {f(d['q1'])} {f(d['q3'])} {f(d['lo'])} {f(d['hi'])}")


# ---------- metrics 2 & 3 (fragment count + MapQ30 rate), both dedup ----------
fc = pd.read_csv(f"{B}/summary/frag_counts_all.tsv", sep="\t")
for col, title in [("uniq_postDedup", "METRIC 2  uniquely-mapped fragments (after dedup)"),
                   ("rate_postDedup", "METRIC 3  MapQ30 rate % (after dedup)"),
                   ("rate_preDedup", "METRIC 3b MapQ30 rate % (before dedup)")]:
    show(title, {t: g[col] for t, g in fc.groupby("dataset")})


# ---------- metrics 4,5,6 (Hi-C contacts) ----------
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

hic = {}   # tech -> list of (cis, cis1kb, trans)
# nagano
rows = []
for f in glob.glob(f"{B}/nagano/alignment/*.summary.txt"):
    d = parse_summary(f)
    rows.append((d.get("UniqMappedMapQ30NoPcrCis"), d.get("UniqMappedMapQ30NoPcrCisMore1kb"),
                 d.get("UniqMappedMapQ30NoPcrTrans")))
hic["nagano"] = rows
# scnomehic (187 passed)
passed = set(pd.read_csv("/home/jmj7858/epifluidlab/scnomehic_paper/s1/QC/gm_passed.txt", header=None)[0].astype(str))
rows = []
for c in passed:
    f = f"{GM}/{c}.summary.txt.gz"
    if os.path.exists(f):
        d = parse_summary(f, gzip.open)
        rows.append((d.get("UniqMappedMapQ30NoPcrCis"), d.get("UniqMappedMapQ30NoPcrCisMore1kb"),
                     d.get("UniqMappedMapQ30NoPcrTrans")))
hic["scnomehic"] = rows
# droplethic (valid.tsv)
dt = pd.read_csv(f"{B}/droplethic/my_project/SRR27586278_hg38.per_cell_qc.valid.tsv", sep="\t")
cc = {c.lower(): c for c in dt.columns}
def col(*names):
    for n in names:
        if n.lower() in cc: return cc[n.lower()]
    return None
c_cis = col("UniqMappedMapQ30NoPcrCis", "UniqMappedMapQ30NoDup_Cis")
c_1kb = col("UniqMappedMapQ30NoPcrCisMore1kb", "UniqMappedMapQ30NoDup_Cis1kb")
c_tr = col("UniqMappedMapQ30NoPcrTrans", "UniqMappedMapQ30NoDup_Trans")
hic["droplethic"] = list(zip(dt[c_cis], dt[c_1kb], dt[c_tr])) if c_cis else []
# snmCseq3 MapQ10 (MappingSummary) and MapQ30 (counts.txt)
ms = pd.read_csv(f"{B}/snmCseq3/alignment/stats/MappingSummary.csv.gz", compression="gzip", index_col=0)
rows10 = []
for _, r in ms.iterrows():
    short, lng, tr = r.get("CisShortContact"), r.get("CisLongContact"), r.get("TransContact")
    if pd.notna(short):
        rows10.append((short + lng, lng, tr))
hic["snmCseq3(MapQ10)"] = rows10
rows30 = []
for f in glob.glob(f"{B}/summary/frag_counts/snmCseq3_contacts_q30/*.q30.counts.txt"):
    d = dict(l.strip().split(",") for l in open(f) if "," in l)
    short, lng, tr = float(d["CisShortContact"]), float(d["CisLongContact"]), float(d["TransContact"])
    rows30.append((short + lng, lng, tr))
hic["snmCseq3(MapQ30)"] = rows30

show("METRIC 4  per-cell cis-contacts", {t: [x[0] for x in v if x[0] is not None] for t, v in hic.items()})
show("METRIC 6  per-cell cis>1kb contacts", {t: [x[1] for x in v if x[1] is not None] for t, v in hic.items()})
show("METRIC 5  trans/cis ratio", {t: [x[2] / (x[0] + x[2]) for x in v if x[0] and x[2] is not None and (x[0] + x[2]) > 0] for t, v in hic.items()})


# ---------- metric 1 (conversion) ----------
cv = pd.read_csv(f"{B}/summary/conversion_percell.csv")
show("METRIC 1  non-CpG conversion %% (chrM)", {t: g["noncpg_chrM"] for t, g in cv.groupby("dataset")})
show("METRIC 1  non-CpG conversion %% (autosome)", {t: g["noncpg_auto"] for t, g in cv.groupby("dataset")})


# ---------- loci ----------
am = pd.read_csv(f"{B}/summary/gch_hcg_counts/all_methods.summary.txt", sep="\t")
show("HCG loci per cell", {t: g["HCG_n"] for t, g in am.groupby("dataset")})
show("GCH loci per cell", {t: g["GCH_n"] for t, g in am.groupby("dataset") if g["GCH_n"].notna().any()})
