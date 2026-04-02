#!/usr/bin/env python3
"""
Sanity check: on bhmem (or any BAM), compare bisulfite-corrected mismatch counts using
  (1) pysam get_aligned_pairs(matches_only=True)
  (2) explicit CIGAR walk (M/=/X only; soft clip and indels excluded from substitution path)

They should agree for typical Illumina CIGARs. Rare ops (P/B) fall back to (1) inside cigar_mx.
"""

from __future__ import annotations

import argparse
import os
import sys

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import (
    count_mismatches_corrected_aligned_pairs,
    count_mismatches_corrected_cigar_mx,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_fasta")
    ap.add_argument("bam", help="e.g. bhmem BAM")
    ap.add_argument("--max-reads", type=int, default=50_000, help="Primary reads to scan")
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    n = 0
    disagree = 0
    max_diff_raw = (0, None)
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        a = count_mismatches_corrected_aligned_pairs(read, fa)
        c = count_mismatches_corrected_cigar_mx(read, fa)
        if a != c:
            disagree += 1
            d = abs(a[0] - c[0])
            if d > max_diff_raw[0]:
                max_diff_raw = (d, read.query_name)
        n += 1
        if args.max_reads and n >= args.max_reads:
            break

    bam.close()
    fa.close()

    print(f"primary_reads_scanned\t{n}")
    print(f"reads_where_methods_differ\t{disagree}")
    print(f"max_abs_diff_raw_between_methods\t{max_diff_raw[0]}\t{max_diff_raw[1] or ''}")
    print(
        "# If counts are identical, aligned_pairs and cigar_mx agree on substitution path."
    )


if __name__ == "__main__":
    main()
