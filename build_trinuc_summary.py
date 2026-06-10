#!/usr/bin/env python3
"""Build a unified trinuc-conversion QC summary across both benchmark arms,
using the SAME pooled context definition (HCG/GCH/noncpg) from the shared
BisSNP parser (parse_bissnp_trinuc.pool). This supersedes the inconsistent
literal-ACG (smallwood_qc.py) and interim Bismark-native (scNOMe) numbers.

Output: long format, one row per (arm, cell, chrom):
    arm, cell_type, CellID, chrom, chrom_role, HCG, HCG_n, GCH, GCH_n, noncpg, noncpg_n

chrM is the cross-arm-comparable conversion control (common to both species).
The autosome (mouse chr19 / human chr21) is arm-specific; HCG there reflects
cell-type biology (ESC vs GM/K562), not just data quality.
"""
import csv
import os
from parse_bissnp_trinuc import pool

BASE = os.path.dirname(os.path.abspath(__file__))
SUF = ".rmdup.RG.trinuc_methy"
CLASSES = ("HCG", "GCH", "noncpg")


def smallwood_cells():
    meta = os.path.join(BASE, "smallwood", "metadata_esc.tsv")
    cells = []
    with open(meta) as fh:
        next(fh)
        for line in fh:
            srr, cond, label = line.rstrip("\n").split("\t")[:3]
            cells.append((srr, cond))   # cond = 2i / Ser
    return cells


def scnome_cells():
    gm = [(f"SRR37296{n}", "GM12878") for n in range(42, 54)]   # SRR3729642..653
    k562 = [(f"K562_{i:02d}", "K562") for i in range(1, 12)]      # K562_01..11
    return gm + k562


ARMS = [
    dict(arm="smallwood_scBS", align_dir=os.path.join(BASE, "smallwood", "05.align_mm10"),
         chroms=[("chrM", "control"), ("chr19", "autosome")], mates=[""],
         cells=smallwood_cells()),
    dict(arm="scNOMe", align_dir=os.path.join(BASE, "scnome", "04.alignment"),
         chroms=[("chrM", "control"), ("chr21", "autosome")], mates=["1", "2"],
         cells=scnome_cells()),
]


def main():
    out_path = os.path.join(BASE, "trinuc_qc_summary.csv")
    rows = []
    for arm in ARMS:
        for cell, ctype in arm["cells"]:
            for chrom, role in arm["chroms"]:
                files = [os.path.join(arm["align_dir"],
                                      f"{cell}{'_'+m if m else ''}{SUF}.{chrom}.txt")
                         for m in arm["mates"]]
                r = pool(files)
                row = dict(arm=arm["arm"], cell_type=ctype, CellID=cell,
                           chrom=chrom, chrom_role=role)
                for cls in CLASSES:
                    n, pct = r[cls]
                    row[cls] = "" if pct is None else round(pct, 4)
                    row[f"{cls}_n"] = "" if n is None else n
                rows.append(row)

    cols = ["arm", "cell_type", "CellID", "chrom", "chrom_role",
            "HCG", "HCG_n", "GCH", "GCH_n", "noncpg", "noncpg_n"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    n_cells = sum(len(a["cells"]) for a in ARMS)
    miss = sum(1 for r in rows if r["HCG"] == "")
    print(f"Wrote {out_path}: {len(rows)} rows ({n_cells} cells x chroms); "
          f"{miss} rows with missing trinuc files")


if __name__ == "__main__":
    main()
