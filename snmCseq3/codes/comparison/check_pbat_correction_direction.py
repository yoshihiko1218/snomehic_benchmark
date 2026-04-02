#!/usr/bin/env python3
"""Check whether bhmem PBAT reads need standard or swapped BS correction.

For each read, compute bisulfite-aware NM two ways:
  - Standard:  forward=skip C(ref)->T(read), reverse=skip G(ref)->A(read)
  - Swapped:   forward=skip G(ref)->A(read), reverse=skip C(ref)->T(read)

Compare both to the bhmem NM:i tag to see which direction is correct.
"""

import sys
import pysam
import numpy as np

BHMEM_BAM = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3/04.bhmem_bam/SRR21549292.bhmem.bam"
GENOME = "/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa"
MAX_READS = 5000


def bs_nm(r, genome_fa, swap=False):
    """Compute bisulfite-aware NM.
    swap=False: standard (fwd: skip C>T, rev: skip G>A)
    swap=True:  PBAT-swapped (fwd: skip G>A, rev: skip C>T)
    """
    seq = r.query_sequence.upper()
    pairs = r.get_aligned_pairs(matches_only=True)
    ref_seq = genome_fa.fetch(
        r.reference_name, r.reference_start, r.reference_end
    ).upper()
    nm = 0
    for qpos, rpos in pairs:
        rb = ref_seq[rpos - r.reference_start]
        qb = seq[qpos]
        if qb == rb:
            continue
        if not swap:
            # standard
            if not r.is_reverse and rb == "C" and qb == "T":
                continue
            if r.is_reverse and rb == "G" and qb == "A":
                continue
        else:
            # swapped for PBAT
            if not r.is_reverse and rb == "G" and qb == "A":
                continue
            if r.is_reverse and rb == "C" and qb == "T":
                continue
        nm += 1
    nm += sum(l for op, l in r.cigartuples if op in (1, 2))
    return nm


def main():
    genome_fa = pysam.FastaFile(GENOME)
    bam = pysam.AlignmentFile(BHMEM_BAM, "rb")

    results = []
    n = 0
    for r in bam:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if not r.has_tag("NM"):
            continue
        tag_nm = r.get_tag("NM")
        is_r2 = bool(r.is_paired and r.is_read2)

        nm_std = bs_nm(r, genome_fa, swap=False)
        nm_swap = bs_nm(r, genome_fa, swap=True)

        results.append({
            "read": r.query_name,
            "is_r2": is_r2,
            "is_reverse": r.is_reverse,
            "tag_nm": tag_nm,
            "nm_std": nm_std,
            "nm_swap": nm_swap,
        })
        n += 1
        if n >= MAX_READS:
            break

    bam.close()
    genome_fa.close()

    # Summarize
    for label, filt in [("R1", False), ("R2", True), ("All", None)]:
        sub = [d for d in results if (filt is None or d["is_r2"] == filt)]
        if not sub:
            continue
        n = len(sub)
        std_eq = sum(1 for d in sub if d["nm_std"] == d["tag_nm"])
        swap_eq = sum(1 for d in sub if d["nm_swap"] == d["tag_nm"])
        std_diff = [abs(d["nm_std"] - d["tag_nm"]) for d in sub]
        swap_diff = [abs(d["nm_swap"] - d["tag_nm"]) for d in sub]
        print("=== %s (n=%d) ===" % (label, n))
        print("  Standard:  exact match to NM:i = %d (%.1f%%),  mean |diff| = %.3f" % (
            std_eq, 100.0 * std_eq / n, np.mean(std_diff)))
        print("  Swapped:   exact match to NM:i = %d (%.1f%%),  mean |diff| = %.3f" % (
            swap_eq, 100.0 * swap_eq / n, np.mean(swap_diff)))
        print()

    # Show first 20 where they disagree
    print("=== First 20 reads where standard != swapped ===")
    print("read\tis_r2\tis_rev\tNM:i\tstd\tswap")
    count = 0
    for d in results:
        if d["nm_std"] != d["nm_swap"]:
            print("%s\t%s\t%s\t%d\t%d\t%d" % (
                d["read"], d["is_r2"], d["is_reverse"],
                d["tag_nm"], d["nm_std"], d["nm_swap"]))
            count += 1
            if count >= 20:
                break


if __name__ == "__main__":
    main()
