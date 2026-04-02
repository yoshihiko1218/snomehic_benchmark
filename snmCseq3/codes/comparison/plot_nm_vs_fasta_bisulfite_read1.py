#!/usr/bin/env python3
"""
Explain why **high fraction NM == recompute** can still yield **moderate Pearson** but **high Spearman**.

Plots read1: ``NM:i`` vs FASTA + bisulfite strand recompute (``count_nm_style_edit_distance``).

Usage:
  python plot_nm_vs_fasta_bisulfite_read1.py ref.fa alignments.bam -o out.pdf [--max-reads 5000]
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
    print("ERROR: pip install matplotlib", file=sys.stderr)
    sys.exit(1)

try:
    from scipy.stats import spearmanr
except ImportError:
    spearmanr = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import count_nm_style_edit_distance


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("-o", "--output", required=True, help="Output PDF path")
    ap.add_argument("--max-reads", type=int, default=5000)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    nm_l, y_l = [], []
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.is_read1:
                continue
            if not read.has_tag("NM"):
                continue
            y = count_nm_style_edit_distance(
                read, fa, True, bisulfite_read2_mode="strand"
            )
            if y < 0:
                continue
            nm_l.append(int(read.get_tag("NM")))
            y_l.append(y)
            if len(nm_l) >= args.max_reads:
                break
    fa.close()

    nm = np.array(nm_l, dtype=float)
    y = np.array(y_l, dtype=float)
    eq = nm == y
    pr = float(np.corrcoef(nm, y)[0, 1])
    sp = float(spearmanr(nm, y).correlation) if spearmanr else float("nan")

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 9))

    ax = axes[0, 0]
    ax.scatter(nm[eq], y[eq], s=10, alpha=0.25, c="#2c5282", label="equal", rasterized=True)
    ax.scatter(nm[~eq], y[~eq], s=14, alpha=0.6, c="#c53030", label="differ", rasterized=True)
    mx = max(nm.max(), y.max(), 1)
    ax.plot([0, mx], [0, mx], "k--", lw=1, alpha=0.45)
    ax.set_xlabel("NM:i")
    ax.set_ylabel("FASTA + bisulfite strand recompute")
    ax.set_title(
        f"Read1 (n={len(nm)})\n"
        f"frac(NM==y)={eq.mean():.3f}  Pearson r={pr:.3f}  Spearman ρ={sp:.3f}"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim(left=-0.5)
    ax.set_ylim(bottom=-0.5)

    ax = axes[0, 1]
    cap = min(50, mx + 2)
    m = (nm <= cap) & (y <= cap)
    ax.scatter(nm[m & eq], y[m & eq], s=8, alpha=0.2, c="#2c5282", rasterized=True)
    ax.scatter(nm[m & ~eq], y[m & ~eq], s=12, alpha=0.55, c="#c53030", rasterized=True)
    ax.plot([0, cap], [0, cap], "k--", lw=1, alpha=0.45)
    ax.set_xlim(-0.5, cap + 0.5)
    ax.set_ylim(-0.5, cap + 0.5)
    ax.set_xlabel("NM:i")
    ax.set_ylabel("Recompute")
    ax.set_title(f"Zoom 0–{int(cap)} (dense region)")

    ax = axes[1, 0]
    d = y - nm
    ax.hist(d, bins=np.arange(d.min() - 0.5, d.max() + 1.5), color="#4a5568", edgecolor="white", lw=0.3)
    ax.axvline(0, color="k", ls="--", alpha=0.5)
    ax.set_xlabel("recompute − NM")
    ax.set_ylabel("count")
    ax.set_title("Residual (most mass at 0; tail hurts Pearson)")

    ax = axes[1, 1]
    ax.axis("off")
    txt = (
        "Why ~90% exact but Pearson ~0.65–0.70?\n\n"
        "• Pearson measures **linear** fit across **all** points.\n"
        "• Most reads sit on y = x (often at low NM), but the **~10%**\n"
        "  that **differ** are **not** scattered along a single line\n"
        "  (recompute can sit **above** NM with variable gap).\n"
        "• Those points add **variance** in (y − NM) without preserving\n"
        "  a tight slope, so **r drops** even though **match rate** is high.\n\n"
        "• **Spearman** (rank) often stays **high (~0.9+)** because ordering\n"
        "  by NM still agrees with ordering by recompute for most reads.\n\n"
        "• **MD-based** recompute tracks NM with Pearson ≈ 1 on read1\n"
        "  (see recompute_nm_style_from_md)."
    )
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left", fontsize=10, family="sans-serif")

    fig.suptitle("NM vs FASTA+bisulfite(strand) — read1", fontsize=12, y=1.02)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Wrote {args.output}")
    print(f"frac_eq\t{eq.mean():.6f}\tPearson\t{pr:.6f}\tSpearman\t{sp:.6f}")


if __name__ == "__main__":
    main()
