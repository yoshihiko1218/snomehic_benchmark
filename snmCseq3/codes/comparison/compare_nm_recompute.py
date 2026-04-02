#!/usr/bin/env python3
"""
Recompute SAM-style edit distance from CIGAR + read + reference FASTA:
  - no_bs:  all substitution mismatches + I/D/N lengths (like NM, genomic)
  - with_bs: same but C/T and G/A conversion positions not counted in substitution part

Compare both to NM:i on primary alignments. Expect no_bs ≈ NM when bases match aligner.
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
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import count_nm_style_edit_distance


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=50_000)
    ap.add_argument(
        "-o",
        "--per-read-tsv",
        default="",
        help="Optional TSV: nm_tag, recompute_no_bs, recompute_with_bs, diff_no_bs",
    )
    ap.add_argument(
        "--pdf",
        default="",
        help="Write comparison scatter plots (NM vs recompute no_bs / with_bs) to this PDF path.",
    )
    ap.add_argument(
        "--bisulfite-read2-mode",
        choices=("strand", "pbat_read2", "pbat_r2_fwd_ga_rev_ct"),
        default="strand",
        help="Mate-2 bisulfite mask when bisulfite-correcting: strand; pbat_read2 (R2 C/T+G/A); "
        "pbat_r2_fwd_ga_rev_ct (R2 forward G/A+A/G, reverse C/T+T/C). See report_nm_recompute_by_mate.py.",
    )
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    nm_tags = []
    no_bs = []
    with_bs = []
    skipped = 0
    n = 0

    fout = open(args.per_read_tsv, "w", newline="") if args.per_read_tsv else None
    if fout:
        fout.write("qname\tnm_tag\trecompute_no_bs\trecompute_with_bs\tdiff_no_bs_nm\n")

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("NM"):
            continue
        nm = int(read.get_tag("NM"))
        a = count_nm_style_edit_distance(read, fa, bisulfite_correct=False)
        b = count_nm_style_edit_distance(
            read,
            fa,
            bisulfite_correct=True,
            bisulfite_read2_mode=args.bisulfite_read2_mode,
        )
        if a < 0:
            skipped += 1
            continue
        nm_tags.append(nm)
        no_bs.append(a)
        with_bs.append(b)
        if fout:
            fout.write(f"{read.query_name}\t{nm}\t{a}\t{b}\t{a - nm}\n")
        n += 1
        if args.max_reads and n >= args.max_reads:
            break

    bam.close()
    fa.close()
    if fout:
        fout.close()

    nm_tags = np.array(nm_tags, dtype=np.int64)
    no_bs = np.array(no_bs, dtype=np.int64)
    with_bs = np.array(with_bs, dtype=np.int64)

    print(f"primary_reads_used\t{len(nm_tags)}")
    print(f"skipped_unsupported_cigar\t{skipped}")
    if len(nm_tags) == 0:
        return

    exact = int((no_bs == nm_tags).sum())
    print(f"fraction_recompute_no_bs_eq_nm\t{exact / len(nm_tags):.6f}\t({exact}/{len(nm_tags)})")
    print(f"mean_abs_diff_no_bs_minus_nm\t{float(np.mean(np.abs(no_bs - nm_tags))):.6f}")
    print(f"max_abs_diff_no_bs_minus_nm\t{int(np.max(np.abs(no_bs - nm_tags)))}")
    print(f"mean_recompute_no_bs\t{float(np.mean(no_bs)):.6f}")
    print(f"mean_recompute_with_bs\t{float(np.mean(with_bs)):.6f}")
    print(f"mean_nm_tag\t{float(np.mean(nm_tags)):.6f}")
    if len(nm_tags) >= 2:
        print(
            f"pearson_nm_vs_recompute_no_bs\t{float(np.corrcoef(nm_tags, no_bs)[0, 1]):.6f}"
        )
        print(
            f"pearson_nm_vs_recompute_with_bs\t{float(np.corrcoef(nm_tags, with_bs)[0, 1]):.6f}"
        )
    print(
        "# recompute_no_bs = genomic subs mismatches + I/D/N lengths (no C/T mask). "
        "Bhmem NM is from the aligner (bisulfite-aware); it often will NOT equal no_bs."
    )

    if args.pdf:
        if plt is None:
            print("ERROR: matplotlib required for --pdf", file=sys.stderr)
            sys.exit(1)
        pr0 = float(np.corrcoef(nm_tags, no_bs)[0, 1]) if len(nm_tags) >= 2 else float("nan")
        pr1 = float(np.corrcoef(nm_tags, with_bs)[0, 1]) if len(nm_tags) >= 2 else float("nan")
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        for ax, y, title, pr in (
            (
                axes[0],
                no_bs,
                "Genomic edit distance (no C/T mask)\n+ indels",
                pr0,
            ),
            (
                axes[1],
                with_bs,
                "Same + bisulfite mask on substitutions\n(C/T, G/A)",
                pr1,
            ),
        ):
            ax.scatter(nm_tags, y, alpha=0.2, s=10, c="#1e3a5f", edgecolors="none", rasterized=True)
            mx = max(float(nm_tags.max()), float(y.max()))
            ax.plot([0, mx], [0, mx], "k--", alpha=0.35, lw=1)
            ax.set_xlabel("NM:i tag")
            ax.set_ylabel("Recomputed distance")
            ax.set_title(f"{title}\nPearson r = {pr:.4f} (n={len(nm_tags)})")
            ax.set_xlim(left=-0.5)
            ax.set_ylim(bottom=-0.5)
        fig.suptitle("NM tag vs recomputed edit distance", fontsize=11, y=1.02)
        fig.tight_layout()
        fig.savefig(args.pdf, dpi=150)
        plt.close(fig)
        print(f"Wrote {args.pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
