#!/usr/bin/env python3
"""
Base-pair level diagnosis when **read 2** ``NM:i`` disagrees with **MD-based** recompute.

Computes the same decomposition as ``count_nm_style_edit_distance_from_md`` but splits:

  - ``n_sub_md``: aligned columns where ``read != MD`` reference base
  - ``n_indel``: insertion + deletion/skip lengths from CIGAR (same as NM-style)

For reads where ``n_sub_md + n_indel != NM:i``, prints example **(ref_pos, ref_md, read, ref_mm10)**
for substitution rows and the **CIGAR** / **MD** tags.

Typical outcomes:

  - **Indel accounting:** rare edge cases where CIGAR ``N``/``P``/padding interacts with NM.
  - **Substitution rows:** ``get_aligned_pairs(with_seq=True)`` vs aligner’s internal NM (rare;
    often off-by-small-integer).

Also compares **read 1** disagreement rate on the same scan for context.

Usage:
  python investigate_read2_nm_md_residual.py ref.fa alignments.bam [--max-reads 8000] \\
    [--detail-reads 12] [--min-diff 1]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter

import os

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import count_nm_style_edit_distance_from_md


def _indel_bases_from_cigar(read) -> int:
    n = 0
    ct = read.cigartuples
    if not ct:
        return 0
    for op, ln in ct:
        if op == 1:
            n += ln
        elif op in (2, 3):
            n += ln
    return n


def _subs_from_md(read) -> tuple[int, list[tuple[int, str, str]]]:
    """(count, list of (ref_pos, ref_md_upper, read_upper)) for mismatches vs MD ref."""
    if read.is_unmapped or read.query_sequence is None or not read.has_tag("MD"):
        return -1, []
    qs = read.query_sequence
    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except (ValueError, AttributeError):
        return -1, []
    rows: list[tuple[int, str, str]] = []
    for t in pairs:
        if len(t) != 3:
            continue
        q, rpos, refb = t
        if q is None or rpos is None or refb is None:
            continue
        ru, rf = qs[q].upper(), str(refb).upper()
        if ru != rf:
            rows.append((int(rpos), rf, ru))
    return len(rows), rows


def _md_recompute_full(read) -> int:
    return count_nm_style_edit_distance_from_md(read)


def _count_md_substitution_letters(md: str) -> int:
    """Count A/C/G/T/N letters in MD outside ``^`` deletion runs (SAM MD mismatch markers)."""
    n = 0
    in_del = False
    for c in md:
        if c == "^":
            in_del = True
            continue
        if in_del:
            if c.isdigit():
                in_del = False
            continue
        if c in "ACGTN":
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Base-level trace: read2 NM vs MD recompute residuals"
    )
    ap.add_argument("reference_fasta", help="mm10.fa (for ref column in detail output)")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=8000, help="Max read2 primaries to scan")
    ap.add_argument(
        "--detail-reads",
        type=int,
        default=10,
        help="How many disagreeing reads to print base-by-base",
    )
    ap.add_argument(
        "--min-diff",
        type=int,
        default=1,
        help="Only detail reads with |md_total - NM| >= this",
    )
    ap.add_argument("--mate", choices=("r2", "r1", "both"), default="r2")
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)

    def scan_mate_fixed(mate: str, limit: int):
        # abs_diff, qname, nm, md_tot, n_sub, n_indel, diff, n_md_letters
        diffs: list[tuple[int, str, int, int, int, int, int, int]] = []
        hist = Counter()
        n_ok = n_bad = n_skip = 0
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for read in bam:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if mate == "r1" and not read.is_read1:
                    continue
                if mate == "r2" and not read.is_read2:
                    continue
                if not read.has_tag("NM"):
                    n_skip += 1
                    continue
                nm = int(read.get_tag("NM"))
                md_tot = _md_recompute_full(read)
                if md_tot < 0:
                    n_skip += 1
                    continue
                n_sub, _ = _subs_from_md(read)
                if n_sub < 0:
                    n_skip += 1
                    continue
                n_indel = _indel_bases_from_cigar(read)
                if n_sub + n_indel != md_tot:
                    print(
                        f"WARN internal sum mismatch {read.query_name}: "
                        f"n_sub+n_indel={n_sub + n_indel} md_tot={md_tot}",
                        file=sys.stderr,
                    )
                d = md_tot - nm
                if d == 0:
                    n_ok += 1
                else:
                    n_bad += 1
                    hist[d] += 1
                    md_tag = read.get_tag("MD")
                    nlet = _count_md_substitution_letters(md_tag)
                    diffs.append((abs(d), read.query_name, nm, md_tot, n_sub, n_indel, d, nlet))
                if n_ok + n_bad >= limit:
                    break
        return {
            "mate": mate,
            "n_ok": n_ok,
            "n_bad": n_bad,
            "n_skip": n_skip,
            "hist": hist,
            "diffs": diffs,
        }

    mates = ("r1", "r2") if args.mate == "both" else (args.mate,)
    summaries = []
    for m in mates:
        summaries.append(scan_mate_fixed(m, args.max_reads))

    for s in summaries:
        tot = s["n_ok"] + s["n_bad"]
        print(f"\n=== {s['mate']} scanned={tot}  NM==MD_recompute {s['n_ok']}  disagree {s['n_bad']}  skip {s['n_skip']} ===")
        if tot:
            print(f"fraction_exact\t{s['n_ok'] / tot:.6f}")
        if s["hist"]:
            print("histogram (md_recompute - NM):", dict(sorted(s["hist"].items())))
        if s["mate"] == "r2" and s["diffs"]:
            close = sum(1 for t in s["diffs"] if abs(t[2] - t[7]) <= 1)
            print(
                f"among_disagreeing_reads\tNM ~= MD_letter_count (|diff|<=1)\t"
                f"{close}/{len(s['diffs'])}"
            )

    # Detail: worst disagreements for read2 (or read1 if --mate r1)
    target_mate = "r2" if args.mate in ("r2", "both") else "r1"
    s2 = (
        summaries[0]
        if len(summaries) == 1
        else next(x for x in summaries if x["mate"] == target_mate)
    )
    diffs = sorted(s2["diffs"], key=lambda x: -x[0])
    print(f"\n## Base-pair detail for up to {args.detail_reads} reads (|diff|>={args.min_diff})")

    shown = 0
    wanted_ordered = [t[1] for t in diffs if abs(t[6]) >= args.min_diff]

    for qn in wanted_ordered:
        if shown >= args.detail_reads:
            break
        read = None
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for r in bam:
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if target_mate == "r1" and not r.is_read1:
                    continue
                if target_mate == "r2" and not r.is_read2:
                    continue
                if r.query_name == qn:
                    read = r
                    break
        if read is None:
            continue
        nm = int(read.get_tag("NM"))
        md_tot = _md_recompute_full(read)
        n_sub, sub_rows = _subs_from_md(read)
        n_indel = _indel_bases_from_cigar(read)
        if abs(md_tot - nm) < args.min_diff:
            continue
        chrom = read.reference_name
        md = read.get_tag("MD") if read.has_tag("MD") else ""
        n_md_letters = _count_md_substitution_letters(md) if md else 0
        print(
            f"\n--- {qn}  {target_mate}  flag={read.flag}  "
            f"NM={nm}  MD_recompute={md_tot}  (n_sub={n_sub} + n_indel={n_indel})  diff={md_tot - nm} ---"
        )
        print(
            f"MD_substitution_letters\t{n_md_letters}\t"
            f"(A/C/G/T/N in MD string; often ~NM when Bhmem encodes conversion sites as letters)"
        )
        print(f"CIGAR\t{read.cigarstring}")
        print(f"MD\t{md[:200]}{'...' if len(md) > 200 else ''}")
        print(
            "substitution_rows (ref_pos 0-based, ref_md, read, ref_mm10)\t"
            f"count={len(sub_rows)}"
        )
        for rpos, rmd, ru in sub_rows[:40]:
            rfa = fa.fetch(chrom, rpos, rpos + 1).upper()
            print(f"  {chrom}:{rpos}\tmd={rmd}\tread={ru}\tmm10={rfa}")
        if len(sub_rows) > 40:
            print(f"  ... {len(sub_rows) - 40} more sub rows")
        shown += 1

    fa.close()
    print(
        "\n# Read2 residuals here are usually NM >> pysam MD recompute with n_sub==0 but many MD "
        "letters: Bhmem's MD uses mismatch letters (e.g. C) at conversion-related loci while "
        "get_aligned_pairs(with_seq=True) still reports read==MD base — so NM tracks letter-count "
        "semantics and pysam's decoder tracks another. Use aligner-reported NM/MD together with "
        "Bhmem docs/source for exact definitions; do not assume standard BWA MD decoding on R2."
    )


if __name__ == "__main__":
    main()
