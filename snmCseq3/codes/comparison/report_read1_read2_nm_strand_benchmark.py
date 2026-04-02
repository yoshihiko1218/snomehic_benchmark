#!/usr/bin/env python3
"""
Compare **read 1 vs read 2** agreement between ``NM:i`` and (1) strand FASTA recompute, (2) MD recompute.

**Bhmem (this project)** is invoked with ``-nonDirectional`` and ``-pbat`` (see ``snmCseq3/codes/02.alignment.sh``):
the library is PBAT + non-directional; mate 1 and mate 2 are not symmetric for bisulfite conversion vs
genomic reference. The strand-only mask is a good **linear** proxy for ``NM`` on read 1 (Spearman
often high) but can be **negatively** correlated with ``NM`` on read 2 — the ordering of reads by
strand-recompute does not match the ordering by ``NM``, not merely a weaker ``r``.

``MD``-based recompute stays much closer to ``NM`` on both mates; read 2 is slightly worse than read 1.

Requires ``scipy`` for Spearman; otherwise only Pearson is printed.

Usage:
  python report_read1_read2_nm_strand_benchmark.py ref.fa alignments.bam [--max-reads 800]
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    from scipy.stats import spearmanr
except ImportError:
    spearmanr = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import (
    count_nm_style_edit_distance,
    count_nm_style_edit_distance_from_md,
)


def _collect(
    bam_path: str,
    fa: pysam.FastaFile,
    mate: str,
    max_reads: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nm_l, rs_l, rmd_l = [], [], []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        n = 0
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if mate == "r1" and not read.is_read1:
                continue
            if mate == "r2" and not read.is_read2:
                continue
            if not read.has_tag("NM"):
                continue
            rs = count_nm_style_edit_distance(
                read, fa, True, bisulfite_read2_mode="strand"
            )
            rmd = count_nm_style_edit_distance_from_md(read)
            if rs < 0 or rmd < 0:
                continue
            nm_l.append(int(read.get_tag("NM")))
            rs_l.append(rs)
            rmd_l.append(rmd)
            n += 1
            if n >= max_reads:
                break
    return (
        np.array(nm_l, dtype=float),
        np.array(rs_l, dtype=float),
        np.array(rmd_l, dtype=float),
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Pearson/Spearman: NM vs strand FASTA vs MD recompute, read1 vs read2."
    )
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=800)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    print(f"max_reads_per_mate\t{args.max_reads}")
    for mate in ("r1", "r2"):
        nm, rs, rmd = _collect(args.bam, fa, mate, args.max_reads)
        if len(nm) < 3:
            print(f"{mate}: too few reads")
            continue
        pr_s = float(np.corrcoef(nm, rs)[0, 1])
        pr_m = float(np.corrcoef(nm, rmd)[0, 1])
        line = (
            f"{mate}\tn={len(nm)}\t"
            f"Pearson(NM,strand)={pr_s:.4f}\tPearson(NM,MD)={pr_m:.4f}\t"
            f"mean|strand-NM|={float(np.mean(np.abs(rs - nm))):.4f}\t"
            f"mean|MD-NM|={float(np.mean(np.abs(rmd - nm))):.4f}"
        )
        if spearmanr is not None:
            sp_s, _ = spearmanr(nm, rs)
            sp_m, _ = spearmanr(nm, rmd)
            line += f"\tSpearman(NM,strand)={float(sp_s):.4f}\tSpearman(NM,MD)={float(sp_m):.4f}"
        print(line)
    fa.close()
    print(
        "\n# Interpretation: negative strand Spearman on read2 means strand FASTA recompute does not "
        "rank reads like NM — use MD recompute or pbat_read2 masking, not strand-only, for mate 2."
    )


if __name__ == "__main__":
    main()
