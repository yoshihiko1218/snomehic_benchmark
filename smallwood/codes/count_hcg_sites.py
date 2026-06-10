#!/usr/bin/env python3
"""Count detected HCG sites per Smallwood ESC cell for the cross-method benchmark.

HCG site count = number of covered positions (rows) in
    06.methy/hcg/<cell>.NOMe.CpG.cov.gz
i.e. the SAME definition scNOMe uses (nome_qc_sites_trinuc.py::count_cov_sites,
union of (chrom,pos)). Smallwood has one cov per cell (no mate split), so this
is just the unique-position count of that single file.

Also reports the all-CpG site count (rows of 06.methy/<cell>.dedup.bismark.cov.gz)
so you can see HCG vs all-CpG side by side.

Usage:
    python codes/count_hcg_sites.py \
        --meta metadata_esc.tsv \
        --methy 06.methy \
        --out 06.methy/hcg/hcg_site_counts.tsv
"""
import argparse
import gzip
import os


def count_cov_sites(path):
    """Unique (chrom,pos) rows in a Bismark .cov(.gz). Returns int or None."""
    if not path or not os.path.exists(path):
        return None
    seen = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            seen.add((f[0], f[1]))
    return len(seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True, help="metadata_esc.tsv (SRR<TAB>condition<TAB>label)")
    ap.add_argument("--methy", default="06.methy", help="06.methy dir (holds all-CpG cov + hcg/)")
    ap.add_argument("--out", required=True, help="output TSV path")
    a = ap.parse_args()

    rows = []
    with open(a.meta) as fh:
        header = fh.readline()  # SRR  condition  label
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            srr, condition, label = parts[0], parts[1], parts[2]
            allcpg = count_cov_sites(os.path.join(a.methy, f"{srr}.dedup.bismark.cov.gz"))
            hcg = count_cov_sites(os.path.join(a.methy, "hcg", f"{srr}.NOMe.CpG.cov.gz"))
            frac = (hcg / allcpg) if (hcg is not None and allcpg) else ""
            rows.append((srr, condition, label, allcpg, hcg, frac))

    with open(a.out, "w") as out:
        out.write("SRR\tcondition\tlabel\tallCpG_site_count\tHCG_site_count\tHCG_over_allCpG\n")
        for srr, condition, label, allcpg, hcg, frac in rows:
            fr = f"{frac:.4f}" if isinstance(frac, float) else ""
            out.write(f"{srr}\t{condition}\t{label}\t{allcpg}\t{hcg}\t{fr}\n")

    n = sum(1 for r in rows if r[4] is not None)
    print(f"Wrote {a.out}: {n}/{len(rows)} cells with HCG counts")


if __name__ == "__main__":
    main()
