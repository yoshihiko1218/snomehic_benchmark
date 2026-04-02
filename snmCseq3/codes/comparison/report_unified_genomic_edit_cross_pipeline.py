#!/usr/bin/env python3
"""
Cross-pipeline **non-MD** comparable edit distance on the **same mm10.fa** + CIGAR walk.

Uses ``count_cross_pipeline_comparable_edit_distance`` (``bisulfite_corrected_mismatch.py``):
  - **Yap/Bismark** (``XR:Z`` = ``CT`` or ``GA``): mask conversion by ``XR``.
  - **Bhmem** (no usable ``XR``): strand on read 1; read 2 uses ``pbat_r2_fwd_ga_rev_ct`` (forward
    G/A+A/G, reverse C/T+T/C), same fallback as :func:`count_cross_pipeline_comparable_edit_distance`.

**Important**
  - This is **one definition** for **fair yap vs bhmem comparison** on genomic reference.
  - It is **not** expected to match **Bismark ``NM``** (that ``NM`` counts raw genomic C vs T mismatches).
  - On **bhmem**, ``NM`` is bisulfite-aware; this walk **approximates** ``NM`` (often closest on read 1).
  - To **validate** each aligner's ``NM`` tag, use ``count_nm_style_edit_distance_from_md`` instead.

Input: ``mapq_comparison.joined.tsv`` (``base_id``, ``is_r1``, …).

Outputs: ``*.per_read.tsv``, ``*.summary.txt``, ``*.plots.pdf``.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

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
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    plt = None
    PdfPages = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import count_cross_pipeline_comparable_edit_distance


def parse_yap_qname(qname: str):
    parts = qname.split("_")
    if len(parts) < 2:
        return None, None
    base, strand = parts[0], parts[1]
    if strand.startswith("1"):
        is_r1 = True
    elif strand.startswith("2"):
        is_r1 = False
    else:
        return None, None
    return base, is_r1


def collect_bhmem_primary(bam_path: str, keys: set) -> dict:
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not (read.flag & 64) and not (read.flag & 128):
                continue
            is_r1 = bool(read.flag & 64)
            key = (read.query_name, is_r1)
            if key not in keys:
                continue
            out[key] = read
    return out


def collect_yap_best(bam_path: str, keys: set) -> dict:
    by_key = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            base, is_r1 = parse_yap_qname(read.query_name)
            if base is None:
                continue
            key = (base, is_r1)
            if key not in keys:
                continue
            parts = read.query_name.split("_")
            is_prim = "-" not in parts[-1]
            mq = read.mapping_quality
            by_key[key].append((read, mq, is_prim))

    best = {}
    for key, cands in by_key.items():
        prim = [x for x in cands if x[2]]
        pool = prim if prim else cands
        best_read = max(pool, key=lambda x: x[1])[0]
        best[key] = best_read
    return best


def load_joined(path: str, max_rows: int) -> list[dict]:
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(
                {
                    "base_id": row["base_id"],
                    "is_r1": int(row["is_r1"]),
                    "mapq_bhmem": int(row["mapq_bhmem"]),
                    "mapq_yap": int(row["mapq_yap"]),
                }
            )
            if max_rows and len(rows) >= max_rows:
                break
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Unified non-MD genomic edit distance (yap XR + bhmem mate rules)."
    )
    ap.add_argument("joined_tsv")
    ap.add_argument("bhmem_bam")
    ap.add_argument("yap_bam")
    ap.add_argument("reference_fasta", help="Genomic FASTA (same mm10 as analyses)")
    ap.add_argument("-o", "--output-prefix", default="unified_genomic_edit_cross")
    ap.add_argument("--max-rows", type=int, default=50_000, help="0 = all joined rows")
    args = ap.parse_args()

    max_rows = args.max_rows if args.max_rows else 0
    rows = load_joined(args.joined_tsv, max_rows)
    keys = {(r["base_id"], bool(r["is_r1"])) for r in rows}

    print(f"Loaded {len(rows)} joined rows", file=sys.stderr)
    bh = collect_bhmem_primary(args.bhmem_bam, keys)
    yp = collect_yap_best(args.yap_bam, keys)
    print(f"  bhmem keys: {len(bh)}  yap keys: {len(yp)}", file=sys.stderr)

    fa = pysam.FastaFile(args.reference_fasta)
    per_read = f"{args.output_prefix}.per_read.tsv"
    summary_path = f"{args.output_prefix}.summary.txt"

    nmb, nmy = [], []
    ub, uy = [], []
    n_ok = 0

    progress_every = 5000

    with open(per_read, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(
            [
                "base_id",
                "is_r1",
                "mapq_bhmem",
                "mapq_yap",
                "nm_bhmem",
                "unified_edit_bhmem",
                "nm_yap",
                "unified_edit_yap",
            ]
        )
        for rec in rows:
            key = (rec["base_id"], rec["is_r1"])
            br = bh.get(key)
            yr = yp.get(key)
            if br is None or yr is None:
                w.writerow(
                    [
                        rec["base_id"],
                        rec["is_r1"],
                        rec["mapq_bhmem"],
                        rec["mapq_yap"],
                        "",
                        "",
                        "",
                        "",
                    ]
                )
                continue
            if not br.has_tag("NM") or not yr.has_tag("NM"):
                w.writerow(
                    [
                        rec["base_id"],
                        rec["is_r1"],
                        rec["mapq_bhmem"],
                        rec["mapq_yap"],
                        "",
                        "",
                        "",
                        "",
                    ]
                )
                continue
            tmb = int(br.get_tag("NM"))
            tmy = int(yr.get_tag("NM"))
            eb = count_cross_pipeline_comparable_edit_distance(br, fa)
            ey = count_cross_pipeline_comparable_edit_distance(yr, fa)
            if eb < 0 or ey < 0:
                w.writerow(
                    [
                        rec["base_id"],
                        rec["is_r1"],
                        rec["mapq_bhmem"],
                        rec["mapq_yap"],
                        tmb,
                        eb if eb >= 0 else "",
                        tmy,
                        ey if ey >= 0 else "",
                    ]
                )
                continue
            w.writerow(
                [
                    rec["base_id"],
                    rec["is_r1"],
                    rec["mapq_bhmem"],
                    rec["mapq_yap"],
                    tmb,
                    eb,
                    tmy,
                    ey,
                ]
            )
            nmb.append(tmb)
            nmy.append(tmy)
            ub.append(eb)
            uy.append(ey)
            n_ok += 1
            if progress_every and n_ok % progress_every == 0:
                print(f"  processed {n_ok} / {len(rows)} paired rows...", file=sys.stderr)

    fa.close()

    nmb = np.array(nmb, dtype=np.int64)
    nmy = np.array(nmy, dtype=np.int64)
    ub = np.array(ub, dtype=np.int64)
    uy = np.array(uy, dtype=np.int64)

    def _concord(tag: np.ndarray, uni: np.ndarray) -> str:
        if len(tag) == 0:
            return "n\t0\n"
        return (
            f"n\t{len(tag)}\n"
            f"fraction_nm_eq_unified\t{float((tag == uni).mean()):.6f}\n"
            f"mean_abs_nm_minus_unified\t{float(np.mean(np.abs(tag - uni))):.6f}\n"
            f"max_abs_nm_minus_unified\t{int(np.max(np.abs(tag - uni)))}\n"
        )

    with open(summary_path, "w") as sf:
        sf.write(f"n_joined_rows\t{len(rows)}\n")
        sf.write(f"n_rows_ok\t{n_ok}\n\n")
        sf.write("## bhmem (unified vs NM tag)\n")
        sf.write(_concord(nmb, ub))
        sf.write("\n## yap (unified vs NM tag; expect poor — Bismark NM is uncorrected)\n")
        sf.write(_concord(nmy, uy))
        if n_ok >= 2:
            sf.write(
                f"\npearson_unified_bhmem_vs_unified_yap\t{float(np.corrcoef(ub, uy)[0, 1]):.6f}\n"
            )
            sf.write(
                f"mean_unified_yap_minus_bhmem\t{float(np.mean(uy - ub)):.6f}\n"
            )
            sf.write(
                f"median_unified_yap_minus_bhmem\t{float(np.median(uy - ub)):.6f}\n"
            )

    print(f"Wrote {per_read}", file=sys.stderr)
    print(f"Wrote {summary_path}", file=sys.stderr)

    if plt is None or PdfPages is None:
        print("WARNING: matplotlib missing, skip PDF", file=sys.stderr)
        return

    pdf_path = f"{args.output_prefix}.plots.pdf"
    pr = lambda a, b: float(np.corrcoef(a, b)[0, 1]) if len(a) >= 2 else float("nan")

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))

    def scatter(ax, x, y, xlab, ylab, title):
        ax.scatter(x, y, alpha=0.12, s=8, c="#1e3a5f", edgecolors="none", rasterized=True)
        if len(x):
            mx = max(float(x.max()), float(y.max()))
            ax.plot([0, mx], [0, mx], "k--", alpha=0.35, lw=1)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.set_xlim(left=-0.5)
        ax.set_ylim(bottom=-0.5)

    scatter(
        axes[0, 0],
        nmb,
        ub,
        "NM tag (bhmem)",
        "Unified genomic edit (bhmem)",
        f"Bhmem\nPearson r = {pr(nmb, ub):.4f} (n={n_ok})",
    )
    scatter(
        axes[0, 1],
        nmy,
        uy,
        "NM tag (yap)",
        "Unified genomic edit (yap)",
        f"Yap (NM uncorrected)\nPearson r = {pr(nmy, uy):.4f} (n={n_ok})",
    )
    scatter(
        axes[1, 0],
        ub,
        uy,
        "Unified (bhmem)",
        "Unified (yap)",
        f"Same definition, both pipelines\nPearson r = {pr(ub, uy):.4f} (n={n_ok})",
    )
    d = uy - ub
    axes[1, 1].hist(
        d,
        bins=min(60, max(10, int(np.ptp(d)) + 2)),
        color="#285e61",
        edgecolor="black",
        alpha=0.85,
    )
    axes[1, 1].set_xlabel("Unified edit (yap) − unified edit (bhmem)")
    axes[1, 1].set_ylabel("Count")
    axes[1, 1].set_title("Cross-pipeline difference")
    fig.suptitle(
        "Non-MD unified genomic edit distance (XR vs bhmem mate rules)",
        fontsize=11,
        y=1.01,
    )
    fig.tight_layout()

    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
    print(f"Wrote {pdf_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
