#!/usr/bin/env python3
"""
Investigate why recomputed genomic edit distance (no C/T mask) is usually > NM:i on bhmem.

Decomposes each read into:
  - n_sub_ct: C/T or G/A conversion-style substitution mismatches vs genome
  - n_sub_other: other substitution mismatches
  - indel_bases: I + D/N lengths

Hypothesis: bhmem NM is bisulfite-aware and does not treat C/T (and rev G/A) like
naive genomic mismatches, so NM ≈ n_sub_other + indel_bases (approximately).
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

from bisulfite_corrected_mismatch import (
    count_nm_style_edit_distance,
    nm_style_distance_breakdown,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=10_000)
    ap.add_argument(
        "--pdf",
        default="",
        help="Optional PDF: residual vs n_sub_ct; NM vs (other+indel)",
    )
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    nm_l, rec_l, ct_l, oth_l, ind_l = [], [], [], [], []
    n = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("NM"):
            continue
        nm = int(read.get_tag("NM"))
        rec = count_nm_style_edit_distance(read, fa, bisulfite_correct=False)
        br = nm_style_distance_breakdown(read, fa)
        if br[0] < 0:
            continue
        n_ct, n_other, n_ind = br
        assert n_ct + n_other + n_ind == rec
        nm_l.append(nm)
        rec_l.append(rec)
        ct_l.append(n_ct)
        oth_l.append(n_other)
        ind_l.append(n_ind)
        n += 1
        if args.max_reads and n >= args.max_reads:
            break

    bam.close()
    fa.close()

    nm_a = np.array(nm_l, dtype=np.int64)
    rec_a = np.array(rec_l, dtype=np.int64)
    ct_a = np.array(ct_l, dtype=np.int64)
    oth_a = np.array(oth_l, dtype=np.int64)
    ind_a = np.array(ind_l, dtype=np.int64)

    other_plus_indel = oth_a + ind_a
    residual = rec_a - nm_a  # recomputed - NM (usually positive)

    print(f"primary_reads\t{len(nm_a)}")
    print(
        f"fraction_recomputed_gt_nm\t{float((rec_a > nm_a).mean()):.4f}\t"
        f"fraction_eq\t{float((rec_a == nm_a).mean()):.4f}"
    )
    print(f"mean_residual_recomputed_minus_nm\t{float(residual.mean()):.4f}")
    print(f"mean_n_sub_ct\t{float(ct_a.mean()):.4f}")
    print(f"mean_n_sub_other\t{float(oth_a.mean()):.4f}")
    print(f"mean_indel_bases\t{float(ind_a.mean()):.4f}")
    print(f"mean_nm\t{float(nm_a.mean()):.4f}")
    print(f"mean_other_plus_indel\t{float(other_plus_indel.mean()):.4f}")

    if len(nm_a) >= 2:
        print(
            f"pearson_residual_vs_n_sub_ct\t{float(np.corrcoef(residual, ct_a)[0, 1]):.6f}"
        )
        print(
            f"pearson_nm_vs_other_plus_indel\t{float(np.corrcoef(nm_a, other_plus_indel)[0, 1]):.6f}"
        )
        print(
            f"mean_abs_nm_minus_other_indel\t{float(np.mean(np.abs(nm_a - other_plus_indel))):.4f}"
        )
        print(
            f"fraction_nm_eq_other_plus_indel\t{float((nm_a == other_plus_indel).mean()):.4f}"
        )

    print(
        "\n# Interpretation: if pearson_residual_vs_n_sub_ct is high and "
        "NM tracks other+indel, bhmem NM largely ignores conversion-style C/T (and rev G/A) "
        "as edit distance vs your naive genomic walk."
    )

    if args.pdf and plt is not None:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        ax = axes[0]
        ax.scatter(ct_a, residual, alpha=0.2, s=8, c="#1a365d", edgecolors="none", rasterized=True)
        ax.set_xlabel("n_sub_ct (C/T or G/A vs genome)")
        ax.set_ylabel("recomputed_no_mask − NM")
        pr = float(np.corrcoef(ct_a, residual)[0, 1]) if len(ct_a) >= 2 else float("nan")
        ax.set_title(f"Residual vs conversion-like subs\nPearson r = {pr:.4f}")

        ax = axes[1]
        ax.scatter(other_plus_indel, nm_a, alpha=0.2, s=8, c="#744210", edgecolors="none", rasterized=True)
        mx = max(float(other_plus_indel.max()), float(nm_a.max()))
        ax.plot([0, mx], [0, mx], "k--", alpha=0.35, lw=1)
        pr2 = (
            float(np.corrcoef(other_plus_indel, nm_a)[0, 1]) if len(nm_a) >= 2 else float("nan")
        )
        ax.set_xlabel("n_sub_other + indel_bases")
        ax.set_ylabel("NM:i tag")
        ax.set_title(f"NM vs (non-conversion subs + indels)\nPearson r = {pr2:.4f}")
        ax.set_xlim(left=-0.5)
        ax.set_ylim(bottom=-0.5)

        fig.suptitle("Why naive genomic recompute > NM (bhmem)", fontsize=11, y=1.02)
        fig.tight_layout()
        fig.savefig(args.pdf, dpi=150)
        plt.close(fig)
        print(f"Wrote {args.pdf}", file=sys.stderr)
    elif args.pdf:
        print("WARNING: matplotlib missing, skip --pdf", file=sys.stderr)


if __name__ == "__main__":
    main()
