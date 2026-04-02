#!/usr/bin/env python3
"""
Deep dive: read1 alignments where strand-masked recompute > NM:i.

Walks the same CIGAR path as count_nm_style_edit_distance, lists each substitution column
(ref, read, genomic pos, whether bisulfite-skipped, whether it contributes to corrected distance)
and indel totals. Helps explain small positive (recompute - NM) gaps on read 1.
"""

from __future__ import annotations

import argparse
import sys

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, __file__.rsplit("/", 1)[0])

from bisulfite_corrected_mismatch import (
    _bisulfite_skip_substitution,
    count_nm_style_edit_distance,
    nm_style_distance_breakdown,
)

_CIGAR_MATCH_MISMATCH_OPS = frozenset((0, 7, 8))


def walk_detail(read, fasta, bisulfite_correct: bool) -> tuple[list[dict], int, int]:
    """Return (substitution_rows, indel_bases, counted_subs)."""
    chrom = read.reference_name
    qs = read.query_sequence
    rev = read.is_reverse
    read2_symmetric = False
    is_read2 = False
    rows: list[dict] = []
    indel_bases = 0
    counted = 0
    q = 0
    r = read.reference_start
    for op, length in read.cigartuples:
        if op in _CIGAR_MATCH_MISMATCH_OPS:
            for _ in range(length):
                rb = qs[q].upper()
                refb = fasta.fetch(chrom, r, r + 1).upper()
                mismatch = rb != refb
                skipped = False
                counts = False
                if mismatch:
                    if bisulfite_correct:
                        skipped = _bisulfite_skip_substitution(
                            refb, rb, rev, read2_symmetric=read2_symmetric, is_read2=is_read2
                        )
                        counts = not skipped
                    else:
                        counts = True
                    if counts:
                        counted += 1
                    rows.append(
                        {
                            "rpos": r,
                            "qpos": q,
                            "ref": refb,
                            "read": rb,
                            "mismatch": True,
                            "bs_skip": skipped,
                            "counts": counts,
                        }
                    )
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
    return rows, indel_bases, counted


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-scan", type=int, default=25_000, help="Max primary read1 to scan")
    ap.add_argument("--top", type=int, default=12, help="Detailed print for this many reads (largest gap)")
    ap.add_argument(
        "-o",
        "--out",
        default="",
        help="Optional text report path (stdout if omitted)",
    )
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    candidates: list[tuple[int, pysam.AlignedSegment, int, int, int, int, int]] = []
    # (delta, read, nm, with_bs, n_ct, n_other, indel)
    scanned = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.is_paired or not read.is_read1:
            continue
        if not read.has_tag("NM"):
            continue
        nm = int(read.get_tag("NM"))
        wb = count_nm_style_edit_distance(
            read, fa, bisulfite_correct=True, bisulfite_read2_mode="strand"
        )
        if wb < 0:
            continue
        ct, oth, ind = nm_style_distance_breakdown(read, fa)
        if ct < 0:
            continue
        delta = wb - nm
        if delta > 0:
            candidates.append((delta, read, nm, wb, ct, oth, ind))
        scanned += 1
        if scanned >= args.max_scan:
            break

    bam.close()
    fa.close()

    candidates.sort(key=lambda t: -t[0])

    lines: list[str] = []
    lines.append(f"scanned_read1_primary\t{scanned}")
    lines.append(f"read1_with_recompute_gt_nm\t{len(candidates)}")
    if not candidates:
        lines.append("# No read1 with corrected recompute > NM in scan range.")
        text = "\n".join(lines) + "\n"
        if args.out:
            open(args.out, "w").write(text)
        else:
            print(text, end="")
        return

    deltas = [t[0] for t in candidates]
    lines.append(f"delta_min\t{min(deltas)}\tdelta_max\t{max(deltas)}")
    lines.append(f"mean_delta_among_positive\t{sum(deltas)/len(deltas):.4f}")

    # Bucket by delta
    from collections import Counter

    c = Counter(deltas)
    lines.append("delta_histogram_top\t" + ", ".join(f"{k}:{c[k]}" for k in sorted(c.keys())[:15]))

    lines.append("")
    lines.append(
        "## Hypothesis check: among positive-delta reads, does gap come from subs vs indels?"
    )
    subs_only = sum(1 for t in candidates if t[5] > 0 and t[6] == 0)
    indel_any = sum(1 for t in candidates if t[6] > 0)  # t[5]=n_other, t[6]=indel
    lines.append(f"positive_delta_with_n_other_ge_1\t{subs_only}")
    lines.append(f"positive_delta_with_indel_gt_0\t{indel_any}")

    lines.append("")
    lines.append(f"## Top {args.top} reads by (corrected_recompute - NM) — substitution detail")
    lines.append("# Columns: rpos is 0-based start on reference (same as BAM POS for first aligned ref base).")

    fa2 = pysam.FastaFile(args.reference_fasta)
    for i, (delta, read, nm, wb, n_ct, n_other, indel) in enumerate(candidates[: args.top]):
        rows, indel_bases, counted_subs = walk_detail(read, fa2, bisulfite_correct=True)
        assert indel_bases == indel, (indel_bases, indel, read.query_name)
        assert counted_subs == n_other, (counted_subs, n_other, read.query_name)
        lines.append("")
        lines.append(f"--- example {i+1} ---")
        lines.append(f"qname\t{read.query_name}")
        lines.append(f"flag\t{read.flag}\treverse\t{read.is_reverse}")
        lines.append(f"cigar\t{read.cigarstring}")
        lines.append(f"NM:i\t{nm}\tstrand_corrected_recompute\t{wb}\tdelta\t{delta}")
        lines.append(
            f"breakdown\tn_CT_or_GA_masked\t{n_ct}\t"
            f"n_other_subs\t{n_other}\tindel_bases\t{indel}"
        )
        if read.has_tag("MD"):
            lines.append(f"MD:Z\t{read.get_tag('MD')}")
        else:
            lines.append("MD:Z\t.")
        counted_rows = [x for x in rows if x["mismatch"] and x["counts"]]
        skipped_rows = [x for x in rows if x["mismatch"] and x["bs_skip"]]
        lines.append(f"mismatch_columns_total\t{len(rows)}")
        lines.append(f"bisulfite_masked_mismatches\t{len(skipped_rows)}")
        lines.append(f"counted_substitutions\t{len(counted_rows)}")
        lines.append("counted_subs_detail\tref\tread\trpos\tqpos")
        for x in counted_rows:
            lines.append(
                f"\t{x['ref']}\t{x['read']}\t{x['rpos']}\t{x['qpos']}"
            )
        if skipped_rows and delta <= 3:
            lines.append("masked_subs_sample\tref\tread\trpos")
            for x in skipped_rows[:8]:
                lines.append(f"\t{x['ref']}\t{x['read']}\t{x['rpos']}")

    fa2.close()

    text = "\n".join(lines) + "\n"
    if args.out:
        open(args.out, "w").write(text)
        print(f"Wrote {args.out}", file=sys.stderr)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
