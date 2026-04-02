#!/usr/bin/env python3
"""
NM:i vs recomputed edit distance, split by mate (read1 / read2) and bisulfite mask mode.

For PBAT / non-directional bisulfite (e.g. bhmem -pbat -nonDirectional), mate 2 often does not
match ``NM`` when using only mapping-strand C/T vs G/A masking; ``pbat_read2`` mode (strand for
R1, symmetric C/T + G/A for R2) reduces the gap. See ``bisulfite_corrected_mismatch`` docstring.
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
    count_nm_style_edit_distance_from_md,
)


def _mate_label(read: pysam.AlignedSegment) -> str:
    if not read.is_paired:
        return "unpaired"
    if read.is_read1:
        return "read1"
    if read.is_read2:
        return "read2"
    return "paired_unknown"


def _stats(nm: np.ndarray, y: np.ndarray) -> dict:
    n = len(nm)
    if n == 0:
        return {
            "n": 0,
            "frac_eq": float("nan"),
            "mean_abs_diff": float("nan"),
            "pearson": float("nan"),
        }
    eq = float((y == nm).mean())
    mad = float(np.mean(np.abs(y - nm)))
    pr = float(np.corrcoef(nm, y)[0, 1]) if n >= 2 else float("nan")
    return {"n": n, "frac_eq": eq, "mean_abs_diff": mad, "pearson": pr}


def _scatter(ax, nm: np.ndarray, y: np.ndarray, title: str, pr: float) -> None:
    ax.scatter(nm, y, alpha=0.2, s=8, c="#1e3a5f", edgecolors="none", rasterized=True)
    if len(nm):
        mx = max(float(nm.max()), float(y.max()))
        ax.plot([0, mx], [0, mx], "k--", alpha=0.35, lw=1)
    ax.set_xlabel("NM:i")
    ax.set_ylabel("Recomputed")
    ax.set_title(f"{title}\nPearson r = {pr:.4f} (n={len(nm)})")
    ax.set_xlim(left=-0.5)
    ax.set_ylim(bottom=-0.5)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="NM vs recompute edit distance by mate (read1/read2) and bisulfite mask."
    )
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=50_000)
    ap.add_argument(
        "-o",
        "--per-read-tsv",
        default="",
        help="TSV: mate, nm, no_bs, with_bs_strand, with_bs_pbat_read2",
    )
    ap.add_argument(
        "--pdf",
        default="",
        help="PDF: 3×3 panels — FASTA recompute (rows 1–2) + MD-based vs NM (row 3) + note.",
    )
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Per mate: lists
    buckets = {
        "read1": {"nm": [], "no_bs": [], "ws": [], "wp": [], "md": []},
        "read2": {"nm": [], "no_bs": [], "ws": [], "wp": [], "md": []},
        "unpaired": {"nm": [], "no_bs": [], "ws": [], "wp": [], "md": []},
        "paired_unknown": {"nm": [], "no_bs": [], "ws": [], "wp": [], "md": []},
    }
    skipped = 0
    skipped_md = 0
    n_used = 0

    fout = open(args.per_read_tsv, "w", newline="") if args.per_read_tsv else None
    if fout:
        fout.write(
            "qname\tmate\tnm_tag\trecompute_no_bs\trecompute_with_bs_strand\t"
            "recompute_with_bs_pbat_read2\trecompute_from_md\n"
        )

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("NM"):
            continue
        nm = int(read.get_tag("NM"))
        a = count_nm_style_edit_distance(read, fa, bisulfite_correct=False)
        b = count_nm_style_edit_distance(read, fa, bisulfite_correct=True, bisulfite_read2_mode="strand")
        c = count_nm_style_edit_distance(read, fa, bisulfite_correct=True, bisulfite_read2_mode="pbat_read2")
        md = count_nm_style_edit_distance_from_md(read)
        if a < 0:
            skipped += 1
            continue
        if md < 0:
            skipped_md += 1

        mate = _mate_label(read)
        key = mate if mate in buckets else "paired_unknown"
        buckets[key]["nm"].append(nm)
        buckets[key]["no_bs"].append(a)
        buckets[key]["ws"].append(b)
        buckets[key]["wp"].append(c)
        buckets[key]["md"].append(md)

        if fout:
            fout.write(f"{read.query_name}\t{mate}\t{nm}\t{a}\t{b}\t{c}\t{md}\n")

        n_used += 1
        if args.max_reads and n_used >= args.max_reads:
            break

    bam.close()
    fa.close()
    if fout:
        fout.close()

    print(f"primary_reads_used\t{n_used}")
    print(f"skipped_unsupported_cigar\t{skipped}")
    print(f"reads_with_md_recompute_failed\t{skipped_md}")

    for name in ("read1", "read2", "unpaired", "paired_unknown"):
        b = buckets[name]
        nm = np.array(b["nm"], dtype=np.int64)
        no_bs = np.array(b["no_bs"], dtype=np.int64)
        ws = np.array(b["ws"], dtype=np.int64)
        wp = np.array(b["wp"], dtype=np.int64)
        md_arr = np.array(b["md"], dtype=np.int64)
        if len(nm) == 0:
            continue
        s0 = _stats(nm, no_bs)
        s1 = _stats(nm, ws)
        s2 = _stats(nm, wp)
        ok = md_arr >= 0
        if ok.any():
            s3 = _stats(nm[ok], md_arr[ok])
        else:
            s3 = {"frac_eq": float("nan"), "mean_abs_diff": float("nan"), "pearson": float("nan")}
        print(f"\n## {name}\tn={s0['n']}")
        print(
            "metric\tfrac_eq_nm\tmean_abs_diff\tpearson_nm_vs_y\t"
            "label"
        )
        print(
            f"recompute_no_bs\t{s0['frac_eq']:.6f}\t{s0['mean_abs_diff']:.6f}\t{s0['pearson']:.6f}\tgenomic subs + indels (FASTA)"
        )
        print(
            f"with_bs_strand\t{s1['frac_eq']:.6f}\t{s1['mean_abs_diff']:.6f}\t{s1['pearson']:.6f}\tC/T or G/A by mapping strand (FASTA)"
        )
        print(
            f"with_bs_pbat_read2\t{s2['frac_eq']:.6f}\t{s2['mean_abs_diff']:.6f}\t{s2['pearson']:.6f}\tR1 strand; R2 symmetric C/T+G/A (FASTA)"
        )
        print(
            f"from_md\t{s3['frac_eq']:.6f}\t{s3['mean_abs_diff']:.6f}\t{s3['pearson']:.6f}\tsubs vs MD ref + indels (matches NM)"
        )

    print(
        "\n# FASTA walks count mismatches vs raw genome; bhmem NM/MD use bisulfite-consistent reference "
        "bases, so use ``from_md`` to validate NM. Strand / pbat_read2 masks are for genomic "
        "interpretation, not bitwise NM reproduction."
    )

    if args.pdf:
        if plt is None:
            print("ERROR: matplotlib required for --pdf", file=sys.stderr)
            sys.exit(1)
        r1 = buckets["read1"]
        r2 = buckets["read2"]
        nm1 = np.array(r1["nm"], dtype=np.int64)
        n1_no = np.array(r1["no_bs"], dtype=np.int64)
        n1_ws = np.array(r1["ws"], dtype=np.int64)
        n1_wp = np.array(r1["wp"], dtype=np.int64)
        n1_md = np.array(r1["md"], dtype=np.int64)
        nm2 = np.array(r2["nm"], dtype=np.int64)
        n2_no = np.array(r2["no_bs"], dtype=np.int64)
        n2_ws = np.array(r2["ws"], dtype=np.int64)
        n2_wp = np.array(r2["wp"], dtype=np.int64)
        n2_md = np.array(r2["md"], dtype=np.int64)

        fig, axes = plt.subplots(3, 3, figsize=(12, 10))
        pr = lambda nm, y: float(np.corrcoef(nm, y)[0, 1]) if len(nm) >= 2 else float("nan")

        _scatter(axes[0, 0], nm1, n1_no, "Read 1: genomic (no mask)", pr(nm1, n1_no))
        _scatter(axes[0, 1], nm1, n1_ws, "Read 1: bisulfite strand mask", pr(nm1, n1_ws))
        _scatter(axes[0, 2], nm1, n1_wp, "Read 1: pbat_read2 mode\n(same as strand for R1)", pr(nm1, n1_wp))
        _scatter(axes[1, 0], nm2, n2_no, "Read 2: genomic (no mask)", pr(nm2, n2_no))
        _scatter(axes[1, 1], nm2, n2_ws, "Read 2: bisulfite strand mask", pr(nm2, n2_ws))
        _scatter(axes[1, 2], nm2, n2_wp, "Read 2: pbat_read2 (R2 symmetric)", pr(nm2, n2_wp))
        ok1 = n1_md >= 0
        ok2 = n2_md >= 0
        _scatter(
            axes[2, 0],
            nm1[ok1],
            n1_md[ok1],
            "Read 1: MD-based recompute\n(vs NM, should be on diagonal)",
            pr(nm1[ok1], n1_md[ok1]) if ok1.any() else float("nan"),
        )
        _scatter(
            axes[2, 1],
            nm2[ok2],
            n2_md[ok2],
            "Read 2: MD-based recompute\n(vs NM)",
            pr(nm2[ok2], n2_md[ok2]) if ok2.any() else float("nan"),
        )
        axes[2, 2].axis("off")
        axes[2, 2].text(
            0.05,
            0.5,
            "FASTA walks compare to raw genome;\nNM/MD use aligner reference\n(bisulfite-consistent).\nUse bottom row to validate NM.",
            fontsize=9,
            va="center",
            transform=axes[2, 2].transAxes,
        )
        fig.suptitle("NM:i vs recomputed edit distance by mate", fontsize=12, y=1.01)
        fig.tight_layout()
        fig.savefig(args.pdf, dpi=150)
        plt.close(fig)
        print(f"Wrote {args.pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
