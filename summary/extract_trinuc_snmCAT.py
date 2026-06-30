#!/usr/bin/env python
"""
Compute per-cell non-CpG / CpG / GpC conversion proxies for snmCAT the SAME way as
the other datasets' trinuc tables (noncpg=ACT%, endo=ACG%, exo=GCT%), but sourced
from the YAP allcools allc files (snmCAT has no BisSNP/bhmem trinuc_methy files).

allc format (tab): chrom, pos, strand, context(4-mer), mc, cov, methylated
The trinucleotide = first 3 chars of the 4-mer context (1 upstream + C + 1 downstream).
Rate = 100 * sum(mc) / sum(cov) over positions whose trinuc matches, per chromosome.
Uses tabix (allc files are .tbi-indexed) to pull only chrM and chr21.

Outputs (matching summary/trinuc/<ds>.<region>.txt; cols sample/noncpg/endo/exo):
  summary/trinuc/snmCAT.chrM.txt
  summary/trinuc/snmCAT.chr21.txt
"""
import glob, os, subprocess, sys

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
ALLC_GLOB = os.path.join(BASE, "snmCAT/mapping_brain/Group*/allc/*.allc.tsv.gz")
TRINUC = {"noncpg": "ACT", "endo": "ACG", "exo": "GCT"}


def rates_for_region(allc, region):
    """Return dict noncpg/endo/exo as % (or '' if no coverage) for one chrom."""
    mc = {k: 0 for k in TRINUC}
    cov = {k: 0 for k in TRINUC}
    try:
        p = subprocess.run(["tabix", allc, region], capture_output=True, text=True)
    except FileNotFoundError:
        sys.exit("tabix not found on PATH")
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        tri = f[3][:3]
        for k, ctx in TRINUC.items():
            if tri == ctx:
                mc[k] += int(f[4])
                cov[k] += int(f[5])
    out = {}
    for k in TRINUC:
        out[k] = (100.0 * mc[k] / cov[k]) if cov[k] > 0 else ""
    return out


def main():
    files = sorted(glob.glob(ALLC_GLOB))
    print(f"snmCAT allc files: {len(files)}")
    regions = {"chrM": "chrM", "chr21": "chr21"}
    out = {r: [] for r in regions}
    for i, allc in enumerate(files, 1):
        sample = os.path.basename(allc)[: -len(".allc.tsv.gz")]
        for r, reg in regions.items():
            v = rates_for_region(allc, reg)
            out[r].append((sample, v["noncpg"], v["endo"], v["exo"]))
        if i % 20 == 0:
            print(f"  {i}/{len(files)}")
    for r in regions:
        path = os.path.join(BASE, "summary/trinuc", f"snmCAT.{r}.txt")
        with open(path, "w") as fh:
            fh.write("sample\tnoncpg\tendo\texo\n")
            for s, a, b, c in out[r]:
                a = f"{a:.3f}" if a != "" else ""
                b = f"{b:.3f}" if b != "" else ""
                c = f"{c:.3f}" if c != "" else ""
                fh.write(f"{s}\t{a}\t{b}\t{c}\n")
        print(f"wrote {path} ({len(out[r])} cells)")


if __name__ == "__main__":
    main()
