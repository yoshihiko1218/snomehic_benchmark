#!/usr/bin/env python3
"""
Explain why **strand** bisulfite masking tracks bhmem ``NM`` well on read 1 but poorly on read 2.

On **read 1**, mapping strand (``FLAG & 16``) lines up with which genomic strand is read, so
masking only C/T (forward) or only G/A (reverse) matches how conversion appears vs the reference.

On **read 2** (PBAT / non-directional libraries), **the same strand rule is wrong**: conversion
still appears as **both** ref C/read T **and** ref G/read A in the genome coordinates, but:

- **Forward** read2: strand masks **C/T** only; **G/A** mismatches are still counted → inflates
  recomputed distance vs bisulfite-aware ``NM`` (Pearson(residual, n_GA) is typically positive).
- **Reverse** read2: strand masks **G/A** only; **C/T** mismatches are still counted → same issue
  (Pearson(residual, n_CT) positive).

``bisulfite_read2_mode="pbat_read2"`` masks **both** classes on mate 2, which partially fixes this
but **does not** match ``NM`` like ``MD``-based recompute: residual also correlates strongly with
**non–(C/T or G/A)** substitution mismatches (``n_other``), i.e. the gap is not purely conversion
polarity.

Usage:
  python diagnose_read2_strand_vs_nm.py ref.fa alignments.bam [--max-reads 8000]
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import (
    _CIGAR_MATCH_MISMATCH_OPS,
    count_nm_style_edit_distance,
    count_nm_style_edit_distance_from_md,
)


def _mx_mismatch_counts(read, fa: pysam.FastaFile) -> tuple[int, int, int, int] | None:
    if read.is_unmapped or read.query_sequence is None:
        return None
    ct = read.cigartuples
    if not ct:
        return None
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return None
    chrom = read.reference_name
    qs = read.query_sequence
    n_ct = n_ga = n_other = 0
    indel_bases = 0
    q = 0
    r = read.reference_start
    for op, length in ct:
        if op in _CIGAR_MATCH_MISMATCH_OPS:
            for _ in range(length):
                rb = qs[q].upper()
                refb = fa.fetch(chrom, r, r + 1).upper()
                if rb != refb:
                    if refb == "C" and rb == "T":
                        n_ct += 1
                    elif refb == "G" and rb == "A":
                        n_ga += 1
                    else:
                        n_other += 1
                q += 1
                r += 1
        elif op == 1:
            indel_bases += length
            q += length
        elif op in (2, 3):
            indel_bases += length
            r += length
        elif op == 4:
            q += length
        elif op == 5:
            pass
    return n_ct, n_ga, n_other, indel_bases


def main() -> None:
    ap = argparse.ArgumentParser(description="Read2 strand vs NM diagnostic (bhmem)")
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=8000, help="Max read2 primaries")
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    rows: list[tuple] = []
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.is_read2 or not read.is_paired:
                continue
            if not read.has_tag("NM"):
                continue
            w = _mx_mismatch_counts(read, fa)
            if w is None:
                continue
            n_ct, n_ga, n_other, indel = w
            nm = int(read.get_tag("NM"))
            rs = count_nm_style_edit_distance(
                read, fa, True, bisulfite_read2_mode="strand"
            )
            rp = count_nm_style_edit_distance(
                read, fa, True, bisulfite_read2_mode="pbat_read2"
            )
            rmd = count_nm_style_edit_distance_from_md(read)
            if rmd < 0 or rs < 0:
                continue
            rows.append(
                (read.is_reverse, nm, rs, rp, rmd, n_ct, n_ga, n_other, indel)
            )
            if len(rows) >= args.max_reads:
                break
    fa.close()

    if not rows:
        print("No read2 records", file=sys.stderr)
        sys.exit(1)

    A = np.array(rows, dtype=float)
    fwd = A[A[:, 0] == 0]
    rev = A[A[:, 0] == 1]

    def block(name: str, blk: np.ndarray) -> None:
        if len(blk) == 0:
            print(f"\n=== {name} n=0 ===")
            return
        nm, rs, rp, rmd = blk[:, 1], blk[:, 2], blk[:, 3], blk[:, 4]
        print(f"\n=== {name} n={len(blk)} ===")
        print(
            f"mean NM {nm.mean():.3f}  strand {rs.mean():.3f}  "
            f"pbat_r2 {rp.mean():.3f}  from_md {rmd.mean():.3f}"
        )
        print(
            f"frac strand==NM {(rs == nm).mean():.4f}  pbat==NM {(rp == nm).mean():.4f}  "
            f"md==NM {(rmd == nm).mean():.4f}"
        )
        print(
            f"mean(strand-NM) {(rs - nm).mean():.3f}  "
            f"mean(pbat-NM) {(rp - nm).mean():.3f}  "
            f"mean(md-NM) {(rmd - nm).mean():.3f}"
        )
        ct, ga, oth = blk[:, 5], blk[:, 6], blk[:, 7]
        print(f"mean raw M/=/X: n_CT {ct.mean():.3f}  n_GA {ga.mean():.3f}  n_other {oth.mean():.3f}")
        res = rs - nm
        if len(blk) > 2:
            print(
                f"Pearson(residual_strand, n_ga) {np.corrcoef(res, ga)[0, 1]:.4f}  "
                f"(n_ct) {np.corrcoef(res, ct)[0, 1]:.4f}  "
                f"(n_other) {np.corrcoef(res, oth)[0, 1]:.4f}"
            )

    block("read2 forward (FLAG rev=False)", fwd)
    block("read2 reverse (FLAG rev=True)", rev)

    print(
        "\nInterpretation: forward R2 strand masks C/T only — residual tends to move with "
        "n_GA; reverse R2 masks G/A only — residual tends to move with n_CT. "
        "Strong n_other correlation means not all gap is conversion polarity; use MD recompute "
        "to match NM."
    )


if __name__ == "__main__":
    main()
