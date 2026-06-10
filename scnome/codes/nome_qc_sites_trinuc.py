#!/usr/bin/env python3
"""
Compute NOMe QC metrics for one cell directly from the Bismark outputs, replacing
the old Bis-tools-derived metrics (6plus2.bed site counts, trinuc_methy files):

  detected sites (union of covered positions across mates 1 and 2):
    HCG_site_count : covered CpG positions      (rows of <cell>_<m>.NOMe.CpG.cov.gz)
    GCH_site_count : covered GpC positions       (rows of <cell>_<m>.NOMe.GpC.cov.gz)

  per-trinucleotide methylation rate on chrM and chr21 (bisulfite-conversion /
  NOMe controls), computed read-level from the deduplicated BAM + reference:
    <chrom>_noncpg (ACT) : a C with prev!=G and next!=G  -> conversion control (~0)
    <chrom>_endo   (ACG) : prev A, next G                -> endogenous CpG methylation
    <chrom>_exo    (GCT) : prev G, next T                -> GpC accessibility (NOMe)

Cytosine strand is taken from Bismark's XG tag (CT = C on + strand; GA = C on -
strand, context read on the reverse complement). Methylation state is the
case of the Bismark XM call (upper = methylated). Positions ignored by --ignore
are not re-applied here; this is a conversion QC summary, not the methylation call.

Usage:
  python nome_qc_sites_trinuc.py --cell K562_01 --methy_dir 05.methy \
      --align_dir 04.alignment --ref <genome.fa> --chroms chrM,chr21 \
      [--bam_suffix .rmdup.bam] --out qc_stats/K562_01.nome_qc.tsv
"""
import argparse
import gzip
import os
from collections import defaultdict

import pysam

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def count_cov_sites(cov_paths):
    """Union of (chrom,pos) across the given .cov.gz files. Returns int or None."""
    seen = set()
    found_any = False
    for p in cov_paths:
        if not os.path.exists(p):
            continue
        found_any = True
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if not line.strip():
                    continue
                c = line.split("\t", 2)
                seen.add((c[0], c[1]))
    return len(seen) if found_any else None


def classify_trinuc(tri):
    """tri is the 3-base context centered on the cytosine, read 5'->3' on the C's
    own strand (so tri[1] == 'C'). Return 'ACT' (non-CpG/non-GpC), 'ACG' (CpG),
    'GCT' (GpC), or None for ambiguous GCG / non-ACGT."""
    if len(tri) != 3 or tri[1] not in "Cc":
        return None
    prev, nxt = tri[0].upper(), tri[2].upper()
    if prev not in "ACGT" or nxt not in "ACGT":
        return None
    is_cpg = nxt == "G"
    is_gpc = prev == "G"
    if is_cpg and is_gpc:
        return None              # GCG: ambiguous, drop (as in the NOMe convention)
    if is_cpg:
        return "ACG"             # endogenous CpG
    if is_gpc:
        return "GCT"             # GpC (NOMe accessibility)
    return "ACT"                 # non-CpG, non-GpC: conversion control


def trinuc_rates(bam_path, ref_fa, chroms):
    """Return {chrom: {ctx: (meth, unmeth)}} for ctx in ACT/ACG/GCT."""
    fa = pysam.FastaFile(ref_fa)
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = {ch: {"ACT": [0, 0], "ACG": [0, 0], "GCT": [0, 0]} for ch in chroms}
    reflen = {c: l for c, l in zip(bam.references, bam.lengths)}
    for chrom in chroms:
        if chrom not in reflen:
            continue
        L = reflen[chrom]
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            tags = dict(read.tags)
            xm = tags.get("XM")
            xg = tags.get("XG")          # 'CT' (+ strand C) or 'GA' (- strand C)
            if xm is None or xg is None:
                continue
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                call = xm[qpos]
                if call in (".", None) or call not in "zZxXhHuU":
                    continue
                if rpos < 1 or rpos > L - 2:
                    continue
                triplet = fa.fetch(chrom, rpos - 1, rpos + 2)  # + strand triplet
                if xg == "GA":
                    # cytosine on - strand; context is revcomp of + strand triplet
                    triplet = revcomp(triplet)
                ctx = classify_trinuc(triplet)
                if ctx is None:
                    continue
                if call.isupper():
                    out[chrom][ctx][0] += 1
                else:
                    out[chrom][ctx][1] += 1
    fa.close()
    bam.close()
    return out


def rate(pair):
    m, u = pair
    tot = m + u
    return (100.0 * m / tot) if tot else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--methy_dir", default="05.methy")
    ap.add_argument("--align_dir", default="04.alignment")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--chroms", default="chrM,chr21")
    ap.add_argument("--bam_suffix", default=".rmdup.bam",
                    help="per-mate BAM suffix; file = <align_dir>/<cell>_<m><suffix>")
    ap.add_argument("--mates", default="1,2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    chroms = [c for c in a.chroms.split(",") if c]
    mates = [m for m in a.mates.split(",") if m]

    cpg = [os.path.join(a.methy_dir, f"{a.cell}_{m}.NOMe.CpG.cov.gz") for m in mates]
    gpc = [os.path.join(a.methy_dir, f"{a.cell}_{m}.NOMe.GpC.cov.gz") for m in mates]
    stats = {
        "CellID": a.cell,
        "HCG_site_count": count_cov_sites(cpg),
        "GCH_site_count": count_cov_sites(gpc),
    }

    # trinuc: pool mates per chrom (sum counts)
    pooled = {ch: {"ACT": [0, 0], "ACG": [0, 0], "GCT": [0, 0]} for ch in chroms}
    for m in mates:
        bam = os.path.join(a.align_dir, f"{a.cell}_{m}{a.bam_suffix}")
        if not os.path.exists(bam):
            print(f"  [WARN] missing BAM {bam}")
            continue
        r = trinuc_rates(bam, a.ref, chroms)
        for ch in chroms:
            for ctx in ("ACT", "ACG", "GCT"):
                pooled[ch][ctx][0] += r[ch][ctx][0]
                pooled[ch][ctx][1] += r[ch][ctx][1]

    label = {"ACT": "noncpg", "ACG": "endo", "GCT": "exo"}
    for ch in chroms:
        for ctx in ("ACT", "ACG", "GCT"):
            stats[f"{ch}_{label[ctx]}"] = rate(pooled[ch][ctx])
            stats[f"{ch}_{label[ctx]}_n"] = sum(pooled[ch][ctx])

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    keys = list(stats.keys())
    with open(a.out, "w") as fh:
        fh.write(",".join(keys) + "\n")
        fh.write(",".join("" if stats[k] is None else str(stats[k]) for k in keys) + "\n")
    print(f"[OK] wrote {a.out}")
    for k in keys:
        print(f"    {k} = {stats[k]}")


if __name__ == "__main__":
    main()
