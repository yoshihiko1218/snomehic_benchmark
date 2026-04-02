#!/usr/bin/env python3
"""
Benchmark **non-MD** FASTA+CIGAR edit distance vs ``NM:i`` under alternative bisulfite skip rules.

**Bhmem.java (PBAT, no -snm3c):** both mates use ``pbat``; read1 is sent to BWA after **G→A**
substitution; read2 after **C→T** substitution. Non-directional pairing considers CT–CT, GA–GA,
CT–GA, GA–CT combinations.

Hypotheses are encoded as ``(read1_skip, read2_skip)`` callables ``(ref, read, is_reverse) -> bool``.

Production non-MD path: ``bisulfite_read2_mode="pbat_r2_fwd_ga_rev_ct"`` in
``bisulfite_corrected_mismatch.py`` matches **h10** (strand read1; read2 forward G/A+A/G, reverse
C/T+T/C), validated vs this script on the same BAM+FASTA.

Run:
  python benchmark_bisulfite_nm_hypotheses.py ref.fa alignments.bam [--max-reads 4000]
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Callable

import numpy as np

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    from scipy.stats import spearmanr
except ImportError:
    print("ERROR: pip install scipy", file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Type: skip(ref_base, read_base, is_reverse) -> True to not count mismatch
SkipFn = Callable[[str, str, bool], bool]


def _count_with_skip(
    read, fasta: pysam.FastaFile, skip_r1: SkipFn, skip_r2: SkipFn
) -> int:
    if read.is_unmapped or read.query_sequence is None:
        return 0
    ct = read.cigartuples
    if not ct:
        return 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return -1
    chrom = read.reference_name
    qs = read.query_sequence
    rev = read.is_reverse
    sk = skip_r1 if read.is_read1 else skip_r2
    subs = 0
    try:
        pairs = read.get_aligned_pairs(matches_only=True)
    except (ValueError, AttributeError):
        return -1
    for q, rpos in pairs:
        if q is None or rpos is None:
            continue
        rb = qs[q].upper()
        refb = fasta.fetch(chrom, rpos, rpos + 1).upper()
        if rb != refb and not sk(refb, rb, rev):
            subs += 1
    indel = 0
    for op, length in ct:
        if op == 1:
            indel += length
        elif op in (2, 3):
            indel += length
    return subs + indel


def _strand_skip(refb: str, rb: str, rev: bool) -> bool:
    if not rev:
        return refb == "C" and rb == "T"
    return refb == "G" and rb == "A"


def _pbat_r2_symmetric(refb: str, rb: str, rev: bool) -> bool:
    return (refb == "C" and rb == "T") or (refb == "G" and rb == "A")


def _all_four_pairs(refb: str, rb: str, rev: bool) -> bool:
    return (refb, rb) in (
        ("C", "T"),
        ("T", "C"),
        ("G", "A"),
        ("A", "G"),
    )


def _r1_ga_only(refb: str, rb: str, rev: bool) -> bool:
    """Java PBAT R1: G→A on read before align — skip genomic G vs read A."""
    return refb == "G" and rb == "A"


def _r1_ct_only(refb: str, rb: str, rev: bool) -> bool:
    return refb == "C" and rb == "T"


def _r1_ga_or_ct(refb: str, rb: str, rev: bool) -> bool:
    return (refb == "G" and rb == "A") or (refb == "C" and rb == "T")


def _r2_ct_or_tc(refb: str, rb: str, rev: bool) -> bool:
    """Java PBAT R2: C→T on read — skip C/T and mirror T/C vs genome."""
    return (refb == "C" and rb == "T") or (refb == "T" and rb == "C")


def _r2_ga_or_ag(refb: str, rb: str, rev: bool) -> bool:
    return (refb == "G" and rb == "A") or (refb == "A" and rb == "G")


def _r1_pbat_java(refb: str, rb: str, rev: bool) -> bool:
    """R1 PBAT: in silico G→A; complement-style on reverse strand often C/T."""
    if not rev:
        return refb == "G" and rb == "A"
    return refb == "C" and rb == "T"


def _r2_pbat_java(refb: str, rb: str, rev: bool) -> bool:
    """R2 PBAT: in silico C→T; on reverse complement expect G/A style."""
    if not rev:
        return (refb == "C" and rb == "T") or (refb == "T" and rb == "C")
    return (refb == "G" and rb == "A") or (refb == "A" and rb == "G")


def _r2_fwd_tc_rev_ga(refb: str, rb: str, rev: bool) -> bool:
    if not rev:
        return (refb, rb) in (("C", "T"), ("T", "C"))
    return (refb, rb) in (("G", "A"), ("A", "G"))


def _r2_fwd_ga_rev_ct(refb: str, rb: str, rev: bool) -> bool:
    if not rev:
        return (refb, rb) in (("G", "A"), ("A", "G"))
    return (refb, rb) in (("C", "T"), ("T", "C"))


def _r2_strand_inverted(refb: str, rb: str, rev: bool) -> bool:
    """Mate-2 strand rule with forward/reverse masks swapped."""
    if not rev:
        return refb == "G" and rb == "A"
    return refb == "C" and rb == "T"


def _r2_fwd_four_rev_cttc(refb: str, rb: str, rev: bool) -> bool:
    if not rev:
        return _all_four_pairs(refb, rb, rev)
    return (refb, rb) in (("C", "T"), ("T", "C"))


def _r2_fwd_cttc_rev_all(refb: str, rb: str, rev: bool) -> bool:
    if not rev:
        return (refb, rb) in (("C", "T"), ("T", "C"))
    return _all_four_pairs(refb, rb, rev)


def _r2_h10_or_strand(refb: str, rb: str, rev: bool) -> bool:
    return _r2_fwd_ga_rev_ct(refb, rb, rev) or _strand_skip(refb, rb, rev)


HYPOTHESES: dict[str, tuple[SkipFn, SkipFn]] = {
    "h0_strand_r2_pbat_sym": (_strand_skip, _pbat_r2_symmetric),
    "h1_strand_r2_fourpair": (_strand_skip, _all_four_pairs),
    "h2_r1_ga_r2_ct_tc": (_r1_ga_only, _r2_ct_or_tc),
    "h3_r1_strand_r2_ct_tc": (_strand_skip, _r2_ct_or_tc),
    "h4_r1_ga_ct_r2_fourpair": (_r1_ga_or_ct, _all_four_pairs),
    "h5_r1_strand_r2_ct_tc_ga_ag": (
        _strand_skip,
        lambda refb, rb, rev: _r2_ct_or_tc(refb, rb, rev)
        or _r2_ga_or_ag(refb, rb, rev),
    ),
    "h6_java_pbat_mate_rules": (_r1_pbat_java, _r2_pbat_java),
    "h7_r1_strand_r2_java_r2": (_strand_skip, lambda r, b, rev: _r2_pbat_java(r, b, rev)),
    "h8_fourpair_both": (_all_four_pairs, _all_four_pairs),
    "h9_strand_r2_fwd_tc_rev_ga": (_strand_skip, _r2_fwd_tc_rev_ga),
    "h10_strand_r2_fwd_ga_rev_ct": (_strand_skip, _r2_fwd_ga_rev_ct),
    "h11_strand_r2_strand_inverted": (_strand_skip, _r2_strand_inverted),
    "h12_strand_r2_fwd4_rev_cttc": (_strand_skip, _r2_fwd_four_rev_cttc),
    "h13_strand_r2_fwd_cttc_rev4": (_strand_skip, _r2_fwd_cttc_rev_all),
    "h14_strand_r2_h10_or_strand": (_strand_skip, _r2_h10_or_strand),
}


def eval_hypothesis(
    bam_path: str,
    fa: pysam.FastaFile,
    mate: str,
    max_reads: int,
    skip_r1: SkipFn,
    skip_r2: SkipFn,
) -> tuple[float, float, int]:
    nm_l, y_l = [], []
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
            y = _count_with_skip(read, fa, skip_r1, skip_r2)
            if y < 0:
                continue
            nm_l.append(int(read.get_tag("NM")))
            y_l.append(y)
            n += 1
            if n >= max_reads:
                break
    if len(nm_l) < 3:
        return float("nan"), float("nan"), 0
    nm = np.array(nm_l, dtype=float)
    y = np.array(y_l, dtype=float)
    sp = float(spearmanr(nm, y).correlation)
    frac = float((nm == y).mean())
    return sp, frac, len(nm)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=4000)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    print(
        f"# Bhmem PBAT: R1 in silico G→A, R2 in silico C→T (Bhmem.java). "
        f"max_reads/mate={args.max_reads}\n"
    )
    rows = []
    for name, (s1, s2) in HYPOTHESES.items():
        sp1, fq1, n1 = eval_hypothesis(args.bam, fa, "r1", args.max_reads, s1, s2)
        sp2, fq2, n2 = eval_hypothesis(args.bam, fa, "r2", args.max_reads, s1, s2)
        rows.append((sp1 + sp2, sp1, fq1, sp2, fq2, name))
        print(
            f"{name}\tR1 Spearman={sp1:.4f}\tfrac_eq={fq1:.4f}\t"
            f"R2 Spearman={sp2:.4f}\tfrac_eq={fq2:.4f}"
        )

    fa.close()
    rows.sort(key=lambda x: -x[0])
    best = rows[0]
    print(
        f"\n# Best sum(R1_sp+R2_sp): {best[-1]}  "
        f"(R1 ρ={best[1]:.4f}, R2 ρ={best[3]:.4f})"
    )


if __name__ == "__main__":
    main()
