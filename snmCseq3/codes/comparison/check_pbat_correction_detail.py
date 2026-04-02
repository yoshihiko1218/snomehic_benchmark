#!/usr/bin/env python3
"""Break down PBAT correction accuracy by mate x strand x direction."""

import pysam
import numpy as np

BHMEM_BAM = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3/04.bhmem_bam/SRR21549292.bhmem.bam"
GENOME = "/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa"
MAX_READS = 5000


def bs_nm(r, genome_fa, swap=False):
    seq = r.query_sequence.upper()
    pairs = r.get_aligned_pairs(matches_only=True)
    ref_seq = genome_fa.fetch(
        r.reference_name, r.reference_start, r.reference_end
    ).upper()
    nm = 0
    ct_skipped = 0
    ga_skipped = 0
    for qpos, rpos in pairs:
        rb = ref_seq[rpos - r.reference_start]
        qb = seq[qpos]
        if qb == rb:
            continue
        if not swap:
            if not r.is_reverse and rb == "C" and qb == "T":
                ct_skipped += 1
                continue
            if r.is_reverse and rb == "G" and qb == "A":
                ga_skipped += 1
                continue
        else:
            if not r.is_reverse and rb == "G" and qb == "A":
                ga_skipped += 1
                continue
            if r.is_reverse and rb == "C" and qb == "T":
                ct_skipped += 1
                continue
        nm += 1
    nm += sum(l for op, l in r.cigartuples if op in (1, 2))
    return nm, ct_skipped, ga_skipped


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

        nm_std, ct_std, ga_std = bs_nm(r, genome_fa, swap=False)
        nm_swap, ct_swap, ga_swap = bs_nm(r, genome_fa, swap=True)

        results.append({
            "is_r2": is_r2,
            "is_reverse": r.is_reverse,
            "tag_nm": tag_nm,
            "nm_std": nm_std,
            "nm_swap": nm_swap,
            "ct_std": ct_std,
            "ga_std": ga_std,
        })
        n += 1
        if n >= MAX_READS:
            break

    bam.close()
    genome_fa.close()

    # Break down by mate x strand
    for mate_label, mate_filt in [("R1", False), ("R2", True)]:
        for strand_label, strand_filt in [("Forward", False), ("Reverse", True)]:
            sub = [d for d in results
                   if d["is_r2"] == mate_filt and d["is_reverse"] == strand_filt]
            if not sub:
                print("%s %s: no reads" % (mate_label, strand_label))
                continue
            ns = len(sub)
            std_eq = sum(1 for d in sub if d["nm_std"] == d["tag_nm"])
            swap_eq = sum(1 for d in sub if d["nm_swap"] == d["tag_nm"])
            std_diff = np.mean([abs(d["nm_std"] - d["tag_nm"]) for d in sub])
            swap_diff = np.mean([abs(d["nm_swap"] - d["tag_nm"]) for d in sub])
            ct_mean = np.mean([d["ct_std"] for d in sub])
            ga_mean = np.mean([d["ga_std"] for d in sub])
            print("=== %s %s (n=%d) ===" % (mate_label, strand_label, ns))
            print("  Standard match:  %d (%.1f%%), mean |diff|=%.3f" % (std_eq, 100*std_eq/ns, std_diff))
            print("  Swapped match:   %d (%.1f%%), mean |diff|=%.3f" % (swap_eq, 100*swap_eq/ns, swap_diff))
            print("  Mean C>T skipped (std): %.1f,  Mean G>A skipped (std): %.1f" % (ct_mean, ga_mean))
            print()


if __name__ == "__main__":
    main()
