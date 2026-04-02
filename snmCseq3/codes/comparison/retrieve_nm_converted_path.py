#!/usr/bin/env python3
"""
Retrieve NM by **fixed CIGAR walk** on bisulfite-**converted** reference (no genome-wide realign).

Uses:

1. **Python** — :func:`bisulfite_corrected_mismatch.count_nm_style_edit_distance_converted_explicit`
   (aligned pairs + indels).
2. **Optional C** — ``cigar_nm_walk`` (same NM definition: ``M`` mismatches + ``X`` + ``I`` + ``D``/``N``).

Pick PBAT trial with :func:`bhmem_equivalent_selection.recompute_nm_bhmem_style_single_pbat`, then
score that trial with both backends and compare to ``NM:i``.

Example::

  make -C snmCseq3/codes/comparison/cigar_nm_walk
  conda activate scnomehic
  python snmCseq3/codes/comparison/retrieve_nm_converted_path.py \\
    snmCseq3/04.bhmem_bam/SRR21549289.bhmem.bam \\
    --c-binary snmCseq3/codes/comparison/cigar_nm_walk/cigar_nm_walk \\
    --max-reads 2000
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

import pysam

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bhmem_equivalent_selection import recompute_nm_bhmem_style_single_pbat  # noqa: E402
from bisulfite_corrected_mismatch import (  # noqa: E402
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
)

DEFAULT_BISULFITE = (
    "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"
)


def nm_via_c_binary(
    exe: str,
    cigar: str,
    query_converted: str,
    ref_span: str,
) -> int | None:
    try:
        out = subprocess.run(
            [exe, cigar, query_converted, ref_span],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    try:
        return int(out.stdout.strip())
    except ValueError:
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", help="BAM (e.g. bhmem)")
    ap.add_argument(
        "bisulfite_genome",
        nargs="?",
        default=DEFAULT_BISULFITE,
        help="Bismark Bisulfite_Genome directory",
    )
    ap.add_argument(
        "--c-binary",
        default="",
        help="Path to cigar_nm_walk executable (optional; build with make -C cigar_nm_walk)",
    )
    ap.add_argument("--max-reads", type=int, default=3000)
    args = ap.parse_args()

    ct_fa = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_fa = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for p in (ct_fa, ga_fa, args.bam):
        if not os.path.isfile(p):
            print(f"ERROR: missing {p}", file=sys.stderr)
            sys.exit(1)

    c_exe = args.c_binary or os.path.join(_SCRIPT_DIR, "cigar_nm_walk", "cigar_nm_walk")
    use_c = os.path.isfile(c_exe) and os.access(c_exe, os.X_OK)

    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)
    bam = pysam.AlignmentFile(args.bam, "rb")

    c_skip = c_eq = c_ne = 0
    tag_eq = tag_ne = 0
    n = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None or not read.cigarstring:
            continue
        if not read.has_tag("NM"):
            continue

        _, label = recompute_nm_bhmem_style_single_pbat(read, fa_ct, fa_ga)
        if label is None:
            continue

        fa = fa_ct if label.genome == "CT" else fa_ga
        cname = bisulfite_converted_contig_name(fa, read.reference_name, label.genome)
        if cname is None:
            continue

        from bisulfite_corrected_mismatch import _pbat_converted_query_variants

        qconv = None
        for tag, qc in _pbat_converted_query_variants(read):
            if tag == label.qtag:
                qconv = qc
                break
        if qconv is None:
            continue

        d_py = count_nm_style_edit_distance_converted_explicit(
            read, fa, ref_contig=cname, query_converted=qconv
        )
        if d_py < 0:
            continue

        rs = int(read.reference_start)
        re_ = int(read.reference_end)
        ref_span = fa.fetch(cname, rs, re_).upper()

        n += 1
        nm_tag = int(read.get_tag("NM"))
        if d_py == nm_tag:
            tag_eq += 1
        else:
            tag_ne += 1

        if use_c:
            nm_c = nm_via_c_binary(c_exe, read.cigarstring, qconv.upper(), ref_span)
            if nm_c is None:
                c_skip += 1
            elif nm_c == d_py:
                c_eq += 1
            else:
                c_ne += 1
                if c_ne <= 3:
                    print(
                        "C!=py",
                        read.query_name,
                        "cigar",
                        read.cigarstring,
                        "py",
                        d_py,
                        "c",
                        nm_c,
                        file=sys.stderr,
                    )
        else:
            c_skip += 1

        if n >= args.max_reads:
            break

    bam.close()
    fa_ct.close()
    fa_ga.close()

    print("BAM:", args.bam)
    print("reads_scored", n)
    print("picked_trial_nm_vs_NM_tag  equal", tag_eq, " not_equal", tag_ne)
    if n:
        print("  frac_tag_eq", tag_eq / n)
    print("Python_explicit vs picked", "ok", n)
    if use_c:
        print("C_walk vs Python  equal", c_eq, " not_equal", c_ne, " skip", c_skip)
        if c_eq + c_ne:
            print("  frac_c_eq_py", c_eq / (c_eq + c_ne))
    else:
        print("C binary: not used (build: make -C", os.path.join(_SCRIPT_DIR, "cigar_nm_walk"), ")")


if __name__ == "__main__":
    main()
