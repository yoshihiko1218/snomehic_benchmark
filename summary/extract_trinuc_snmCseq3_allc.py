#!/usr/bin/env python
"""
Recompute snmCseq3 (snm3C-seq) per-cell bisulfite-conversion trinuc proxy from its YAP
allc files, POSITION-LEVEL, to replace the noisy read-level bhmem chrM values (only ~7
ACT positions/cell on the tiny chrM -> mostly 0).

snm3C-seq allc context is a 3-mer [C][down1][down2] with NO upstream base, so we recover
the exact ACT/ACG/GCT trinucleotide by looking up the mm10 reference (strand-aware):
  window = ref[pos-1 .. pos+1] (1-based, + strand); trinuc = window if '+' else revcomp.
  noncpg=ACT%, endo=ACG%, exo=GCT%  (= 100 * sum(mc) / sum(cov) per trinuc)
Writes (overwrites) summary/trinuc/snmCseq3.chrM.txt and snmCseq3.chr21.txt (autosome=chr19).
"""
import glob, os, subprocess
import pysam

B = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
MM10 = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa"
TRI = {"noncpg": "ACT", "endo": "ACG", "exo": "GCT"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rates(allc, chrom, seq):
    # seq = whole-chromosome uppercase string (0-based); slice is O(1) per position.
    mc = {v: 0 for v in TRI.values()}
    cov = {v: 0 for v in TRI.values()}
    p = subprocess.run(["tabix", allc, chrom], capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        pos = int(f[1]); strand = f[2]; m = int(f[4]); c = int(f[5])
        w = seq[pos - 2:pos + 1]   # 1-based pos-1,pos,pos+1
        if len(w) != 3:
            continue
        tri = w if strand == "+" else revcomp(w)
        if tri[1] != "C":
            continue
        if tri in mc:
            mc[tri] += m; cov[tri] += c
    out = {}
    for k, t in TRI.items():
        out[k] = (100.0 * mc[t] / cov[t]) if cov[t] > 0 else ""
    return out


def main():
    ref = pysam.FastaFile(MM10)
    allcs = sorted(glob.glob(f"{B}/snmCseq3/alignment/Group*/allc/*.allc.tsv.gz"))
    print(f"snmCseq3 allc files: {len(allcs)}")
    # region -> output file (autosome chr19 saved as *.chr21.txt for collect_conversion.py)
    regions = {"chrM": f"{B}/summary/trinuc/snmCseq3.chrM.txt",
               "chr19": f"{B}/summary/trinuc/snmCseq3.chr21.txt"}
    seqs = {r: ref.fetch(r).upper() for r in regions}   # load each chromosome once
    print("loaded chrom seqs:", {r: len(s) for r, s in seqs.items()})
    rows = {r: [] for r in regions}
    for i, allc in enumerate(allcs, 1):
        cell = os.path.basename(allc)[: -len(".allc.tsv.gz")]
        for r in regions:
            v = rates(allc, r, seqs[r])
            rows[r].append((cell, v["noncpg"], v["endo"], v["exo"]))
        if i % 20 == 0:
            print(f"  {i}/{len(allcs)}")
    for r, path in regions.items():
        with open(path, "w") as fh:
            fh.write("sample\tnoncpg\tendo\texo\n")
            for s, a, b, c in rows[r]:
                fmt = lambda x: f"{x:.3f}" if x != "" else ""
                fh.write(f"{s}\t{fmt(a)}\t{fmt(b)}\t{fmt(c)}\n")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
