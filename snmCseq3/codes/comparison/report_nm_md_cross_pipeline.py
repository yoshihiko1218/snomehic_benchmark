#!/usr/bin/env python3
"""
Cross-pipeline metric: **one method, two BAMs** — recompute edit distance from ``MD:Z`` + CIGAR indels.

This is the practical answer to “first get bhmem mismatch close to ``NM:i``, then apply the **same**
thing on yap”: call ``recompute_nm_style_from_md`` / ``count_nm_style_edit_distance_from_md`` on
every read (``bisulfite_corrected_mismatch.py``). No FASTA, no pipeline-specific masking branches.

- **Bhmem:** recomputed value typically tracks ``NM:i`` very closely (``MD`` matches the aligner’s
  bisulfite-consistent encoding).
- **Yap/Bismark:** recomputed value typically **equals** ``NM:i`` (``NM`` and ``MD`` are tied).

Pooled read1+read2 on bhmem can show a small residual vs ``NM``; split by mate in
``report_nm_recompute_by_mate.py`` if needed. **Raw mm10 FASTA** walks do not approximate bhmem
``NM`` as well as ``MD``-based recompute.

Input: a **joined** TSV with columns ``base_id``, ``is_r1`` (e.g. ``mapq_comparison.joined.tsv`` from
``compare_mapq.py``). For each row, loads the primary bhmem alignment and the best yap alignment
(same logic as ``report_discrepant_mismatch.py``).

Outputs:
  - ``*.per_read.tsv`` — tags + MD-recompute for both
  - ``*.summary.txt`` — concordance ``NM`` vs MD-recompute per BAM; correlation MD-recomputes
  - ``*.plots.pdf`` — scatter ``NM`` vs MD (each BAM) + cross scatter MD_bhmem vs MD_yap + diff histogram
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

from bisulfite_corrected_mismatch import recompute_nm_style_from_md


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
    ap = argparse.ArgumentParser(description="Cross-pipeline NM via MD-recompute (bhmem + yap).")
    ap.add_argument("joined_tsv", help="mapq_comparison.joined.tsv (base_id, is_r1, mapq_*)")
    ap.add_argument("bhmem_bam")
    ap.add_argument("yap_bam")
    ap.add_argument("-o", "--output-prefix", default="nm_md_cross_pipeline")
    ap.add_argument("--max-rows", type=int, default=50_000, help="Max joined rows (0 = all)")
    args = ap.parse_args()

    max_rows = args.max_rows if args.max_rows else 0
    rows = load_joined(args.joined_tsv, max_rows)
    keys = {(r["base_id"], bool(r["is_r1"])) for r in rows}

    print(f"Loaded {len(rows)} joined rows", file=sys.stderr)
    print("Collecting bhmem...", file=sys.stderr)
    bh = collect_bhmem_primary(args.bhmem_bam, keys)
    print(f"  bhmem keys found: {len(bh)}", file=sys.stderr)
    print("Collecting yap...", file=sys.stderr)
    yp = collect_yap_best(args.yap_bam, keys)
    print(f"  yap keys found: {len(yp)}", file=sys.stderr)

    per_read = f"{args.output_prefix}.per_read.tsv"
    summary_path = f"{args.output_prefix}.summary.txt"

    nmb_tag, nmb_md = [], []
    nmy_tag, nmy_md = [], []
    n_ok = 0

    with open(per_read, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(
            [
                "base_id",
                "is_r1",
                "mapq_bhmem",
                "mapq_yap",
                "nm_bhmem_tag",
                "nm_bhmem_from_md",
                "nm_yap_tag",
                "nm_yap_from_md",
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
            mdb = recompute_nm_style_from_md(br)
            mdy = recompute_nm_style_from_md(yr)
            if mdb < 0 or mdy < 0:
                w.writerow(
                    [
                        rec["base_id"],
                        rec["is_r1"],
                        rec["mapq_bhmem"],
                        rec["mapq_yap"],
                        tmb,
                        mdb if mdb >= 0 else "",
                        tmy,
                        mdy if mdy >= 0 else "",
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
                    mdb,
                    tmy,
                    mdy,
                ]
            )
            nmb_tag.append(tmb)
            nmb_md.append(mdb)
            nmy_tag.append(tmy)
            nmy_md.append(mdy)
            n_ok += 1

    nmb_tag = np.array(nmb_tag, dtype=np.int64)
    nmb_md = np.array(nmb_md, dtype=np.int64)
    nmy_tag = np.array(nmy_tag, dtype=np.int64)
    nmy_md = np.array(nmy_md, dtype=np.int64)

    def _stats(tag: np.ndarray, md: np.ndarray) -> str:
        if len(tag) == 0:
            return "n\t0\n"
        eq = float((tag == md).mean())
        mad = float(np.mean(np.abs(tag - md)))
        mx = int(np.max(np.abs(tag - md)))
        return (
            f"n\t{len(tag)}\n"
            f"fraction_nm_tag_eq_nm_from_md\t{eq:.6f}\n"
            f"mean_abs_nm_tag_minus_nm_from_md\t{mad:.6f}\n"
            f"max_abs_nm_tag_minus_nm_from_md\t{mx}\n"
        )

    with open(summary_path, "w") as sf:
        sf.write(f"n_joined_rows\t{len(rows)}\n")
        sf.write(f"n_rows_with_both_bams_and_nm_and_md\t{n_ok}\n\n")
        sf.write("## bhmem\n")
        sf.write(_stats(nmb_tag, nmb_md))
        sf.write("\n## yap\n")
        sf.write(_stats(nmy_tag, nmy_md))
        if n_ok >= 2:
            pr = float(np.corrcoef(nmb_md, nmy_md)[0, 1])
            sf.write(f"\npearson_nm_from_md_bhmem_vs_yap\t{pr:.6f}\n")
            sf.write(
                f"mean_nm_from_md_yap_minus_bhmem\t{float(np.mean(nmy_md - nmb_md)):.6f}\n"
            )
            sf.write(
                f"median_nm_from_md_yap_minus_bhmem\t{float(np.median(nmy_md - nmb_md)):.6f}\n"
            )

    print(f"Wrote {per_read}", file=sys.stderr)
    print(f"Wrote {summary_path}", file=sys.stderr)

    if plt is None or PdfPages is None:
        print("WARNING: matplotlib missing, skip PDF", file=sys.stderr)
        return

    pdf_path = f"{args.output_prefix}.plots.pdf"
    pr = lambda a, b: float(np.corrcoef(a, b)[0, 1]) if len(a) >= 2 else float("nan")

    fig1, axes = plt.subplots(2, 2, figsize=(10, 9))

    def scatter(ax, x, y, xlab, ylab, title):
        ax.scatter(x, y, alpha=0.15, s=8, c="#1e3a5f", edgecolors="none", rasterized=True)
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
        nmb_tag,
        nmb_md,
        "NM tag (bhmem)",
        "NM from MD (bhmem)",
        f"Bhmem: tag vs MD-recompute\nPearson r = {pr(nmb_tag, nmb_md):.4f} (n={n_ok})",
    )
    scatter(
        axes[0, 1],
        nmy_tag,
        nmy_md,
        "NM tag (yap)",
        "NM from MD (yap)",
        f"Yap: tag vs MD-recompute\nPearson r = {pr(nmy_tag, nmy_md):.4f} (n={n_ok})",
    )
    scatter(
        axes[1, 0],
        nmb_md,
        nmy_md,
        "NM from MD (bhmem)",
        "NM from MD (yap)",
        f"Cross-pipeline (same definition)\nPearson r = {pr(nmb_md, nmy_md):.4f} (n={n_ok})",
    )
    d = nmy_md - nmb_md
    axes[1, 1].hist(d, bins=min(60, max(10, int(np.ptp(d)) + 2)), color="#553c9a", edgecolor="black", alpha=0.85)
    axes[1, 1].set_xlabel("NM from MD (yap) − NM from MD (bhmem)")
    axes[1, 1].set_ylabel("Count")
    axes[1, 1].set_title("Difference (MD-based metric)")
    fig1.suptitle("Consistent NM: MD + CIGAR indels (both pipelines)", fontsize=12, y=1.01)
    fig1.tight_layout()

    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig1, dpi=150)
        plt.close(fig1)
    print(f"Wrote {pdf_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
