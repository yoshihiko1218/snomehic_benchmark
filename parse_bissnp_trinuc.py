#!/usr/bin/env python3
"""Shared parser for BisSNP/Bis-tools `*.trinuc_methy.<chrom>.txt` files.

Each trinuc_methy file reports per-trinucleotide methylation (16 contexts), e.g.
    ACG:    10753   77.606%
This parser pools them into the NOMe context CLASSES, consistently for any arm:
    HCG    (endogenous CpG meth) = ACG + CCG + TCG          (next=G, prev!=G)
    GCH    (GpC accessibility)   = GCA + GCC + GCT          (prev=G, next!=G)
    noncpg (conversion control)  = ACA ACC ACT CCA CCC CCT TCA TCC TCT
    GCG is dropped (ambiguous: both CpG and GpC) — same convention as Bismark
    coverage2cytosine --nome-seq and nome_qc_sites_trinuc.py.

Pooling is by COUNTS: meth_i = n_i * pct_i/100; class% = 100*sum(meth)/sum(n).
This is the SAME context definition as the Bismark-native route, so BisSNP and
Bismark-native numbers are finally on the same footing — the only remaining
differences are BisSNP's SNP/MAPQ/coverage filtering (intended: cleaner calls).

CLI (one cell; pools across mates and writes one CSV row):
    python parse_bissnp_trinuc.py --cell SRR1248481 \
        --align_dir smallwood/05.align_mm10 --chroms chrM,chr19 --mates "" \
        --out qc_stats/bissnp_trinuc/SRR1248481.csv
"""
import argparse
import math
import os

HCG = {"ACG", "CCG", "TCG"}
GCH = {"GCA", "GCC", "GCT"}
NONCPG = {"ACA", "ACC", "ACT", "CCA", "CCC", "CCT", "TCA", "TCC", "TCT"}
# GCG intentionally excluded (ambiguous)


def parse_file(path):
    """Return {trinuc: (n, meth)} for a single trinuc_methy file, or {} if absent."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            tri = parts[0].rstrip(":").upper()
            if len(tri) != 3:
                continue
            try:
                n = int(parts[1])
                pct = float(parts[2].rstrip("%"))
            except ValueError:
                continue
            # BisSNP writes "NaN%" for zero-coverage contexts; skip them
            # (0 * NaN = NaN would otherwise poison the pooled sum).
            if n <= 0 or not math.isfinite(pct):
                continue
            out[tri] = (n, n * pct / 100.0)
    return out


def pool(files):
    """Pool a list of trinuc_methy files (e.g. mates) into class-level (n, %).
    Returns {'HCG': (n, pct|None), 'GCH': (...), 'noncpg': (...)}; pct None if n=0,
    and the whole value None if NO input file existed."""
    classes = {"HCG": HCG, "GCH": GCH, "noncpg": NONCPG}
    tot = {k: [0, 0.0] for k in classes}            # [n, meth]
    found = False
    for p in files:
        d = parse_file(p)
        if d:
            found = True
        for tri, (n, m) in d.items():
            for cname, cset in classes.items():
                if tri in cset:
                    tot[cname][0] += n
                    tot[cname][1] += m
    if not found:
        return {k: (None, None) for k in classes}
    res = {}
    for k, (n, m) in tot.items():
        res[k] = (n, (100.0 * m / n) if n else None)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--align_dir", required=True)
    ap.add_argument("--trinuc_suffix", default=".rmdup.RG.trinuc_methy")
    ap.add_argument("--chroms", default="chrM,chr21")
    ap.add_argument("--mates", default="1,2", help='comma list; "" = single no-mate cell')
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    chroms = [c for c in a.chroms.split(",") if c]
    mates = [m for m in a.mates.split(",") if m] or [""]

    stats = {"CellID": a.cell}
    for ch in chroms:
        files = [os.path.join(a.align_dir,
                              f"{a.cell}{'_'+m if m else ''}{a.trinuc_suffix}.{ch}.txt")
                 for m in mates]
        r = pool(files)
        for cls in ("HCG", "GCH", "noncpg"):
            n, pct = r[cls]
            stats[f"{ch}_{cls}"] = "" if pct is None else f"{pct:.4f}"
            stats[f"{ch}_{cls}_n"] = "" if n is None else n

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    keys = list(stats.keys())
    with open(a.out, "w") as fh:
        fh.write(",".join(keys) + "\n")
        fh.write(",".join(str(stats[k]) for k in keys) + "\n")
    print(f"[OK] {a.cell} -> {a.out}")


if __name__ == "__main__":
    main()
