#!/usr/bin/env python3
"""
For per_read rows with nm_bhmem==0 and yap_corrected_mismatch==0, compare
bhmem vs yap primary/best alignment coordinates (chrom, start, end, strand).

Also reports: optional same-strand locus match within --tolerance bp (max of
start/end absolute differences), and summary stats for max(|Δstart|,|Δend|) on
same chromosome when spans differ.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict

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


def collect_bhmem_loc(bam_path: str, keys: set) -> dict:
    """key -> (ref, start, end, is_reverse) primary alignment."""
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not (read.flag & 64) and not (read.flag & 128):
                continue
            qname = read.query_name
            is_r1 = bool(read.flag & 64)
            key = (qname, is_r1)
            if key not in keys:
                continue
            out[key] = (
                read.reference_name,
                read.reference_start,
                read.reference_end,
                read.is_reverse,
            )
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


def write_report_pdf(
    pdf_path: str,
    *,
    per_read_tsv: str,
    n_total: int,
    same_strict: int,
    same_chrom_diff_span: int,
    diff_chrom: int,
    strand_mismatch: int,
    missing_bh: int,
    missing_yap: int,
    tol: int,
    within_tol: int,
    max_deltas_diff_span: list[int],
) -> None:
    if plt is None or PdfPages is None:
        print("ERROR: matplotlib required for --pdf", file=sys.stderr)
        sys.exit(1)

    labels = [
        "Same locus (strict)",
        "Same chr, diff span",
        "Different chr",
        "Same span, diff strand",
    ]
    counts = [
        same_strict,
        same_chrom_diff_span,
        diff_chrom,
        strand_mismatch,
    ]
    colors = ["#2ecc71", "#f39c12", "#e74c3c", "#9b59b6"]

    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(2, 1, figsize=(8.5, 10))
        fig.suptitle(
            "Bhmem vs Yap alignment (NM_bhmem=0, yap corrected mismatch=0)",
            fontsize=13,
            fontweight="bold",
        )

        ax0 = axes[0]
        y_pos = range(len(labels))
        ax0.barh(y_pos, counts, color=colors, edgecolor="white", linewidth=0.8)
        ax0.set_yticks(list(y_pos))
        ax0.set_yticklabels(labels)
        ax0.set_xlabel("Read count")
        ax0.set_title("Agreement on primary alignment coordinates")
        ax0.invert_yaxis()
        for i, c in enumerate(counts):
            ax0.text(
                c + max(counts) * 0.01,
                i,
                str(c),
                va="center",
                fontsize=10,
            )

        ax1 = axes[1]
        ax1.axis("off")
        lines = [
            f"Source per-read TSV: {per_read_tsv}",
            f"Reads (zero-zero): {n_total}",
            "",
            "Strict match: same chr, start, end, strand.",
            f"Within {tol} bp (same chr + strand, max |Δstart|,|Δend| ≤ {tol}): {within_tol}",
            "",
            f"Missing bhmem: {missing_bh}  |  Missing yap: {missing_yap}",
        ]
        if max_deltas_diff_span:
            lines.extend(
                [
                    "",
                    f"Same chr, different span (n={len(max_deltas_diff_span)}):",
                    f"  max(|Δstart|,|Δend|) min = {min(max_deltas_diff_span)} bp",
                    f"  max(|Δstart|,|Δend|) median = {statistics.median(max_deltas_diff_span):.0f} bp",
                    f"  max(|Δstart|,|Δend|) max = {max(max_deltas_diff_span)} bp",
                    "",
                    "Large medians usually mean different loci on the same chromosome,",
                    "not soft-clip noise.",
                ]
            )
        ax1.text(
            0.02,
            0.98,
            "\n".join(lines),
            transform=ax1.transAxes,
            fontsize=9,
            verticalalignment="top",
            fontfamily="monospace",
        )

        fig.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.06, hspace=0.35)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        if max_deltas_diff_span:
            fig2, ax = plt.subplots(figsize=(8.5, 5))
            logs = [math.log10(x + 1) for x in max_deltas_diff_span]
            ax.hist(logs, bins=40, color="#3498db", edgecolor="white", linewidth=0.5)
            ax.set_xlabel("log10(max(|Δstart|, |Δend|) + 1)")
            ax.set_ylabel("Reads")
            ax.set_title(
                "Same chromosome, different span: size of coordinate difference (yap vs bhmem)"
            )
            plt.tight_layout()
            pdf.savefig(fig2, bbox_inches="tight")
            plt.close(fig2)

            small = [x for x in max_deltas_diff_span if x <= 10_000]
            if len(small) >= 2 and len(small) < len(max_deltas_diff_span):
                fig3, ax = plt.subplots(figsize=(8.5, 5))
                ax.hist(small, bins=min(50, max(10, len(small) // 5)), color="#16a085", edgecolor="white", linewidth=0.5)
                ax.set_xlabel("max(|Δstart|, |Δend|) [bp], ≤ 10,000 only")
                ax.set_ylabel("Reads")
                ax.set_title(
                    f"Local differences only (n={len(small)} / {len(max_deltas_diff_span)} same-chr-diff-span reads)"
                )
                plt.tight_layout()
                pdf.savefig(fig3, bbox_inches="tight")
                plt.close(fig3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("per_read_tsv", help="discrepant_mismatch_report*.per_read.tsv")
    ap.add_argument("bhmem_bam")
    ap.add_argument("yap_bam")
    ap.add_argument(
        "-o",
        "--output-tsv",
        default="",
        help="Optional path for per-read comparison table",
    )
    ap.add_argument(
        "--tolerance",
        type=int,
        default=5,
        help="Same chrom + same strand: count as matching locus if "
        "max(|yap_start-bhmem_start|, |yap_end-bhmem_end|) <= this (bp). Default 5.",
    )
    ap.add_argument(
        "--pdf",
        dest="pdf_path",
        default="",
        help="Write summary figures to this PDF path (requires matplotlib).",
    )
    args = ap.parse_args()
    tol = args.tolerance

    zero_rows = []
    with open(args.per_read_tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                nmb = row.get("nm_bhmem", "").strip()
                yc = row.get("yap_corrected_mismatch", "").strip()
                if nmb == "" or yc == "":
                    continue
                if int(nmb) != 0 or int(yc) != 0:
                    continue
            except ValueError:
                continue
            bid = row["base_id"]
            is_r1 = bool(int(row["is_r1"]))
            zero_rows.append((bid, is_r1))

    keys = set(zero_rows)
    print(f"rows with nm_bhmem=0 and yap_corrected_mismatch=0: {len(zero_rows)}", file=sys.stderr)
    print("Collecting bhmem coordinates...", file=sys.stderr)
    bh_loc = collect_bhmem_loc(args.bhmem_bam, keys)
    print("Collecting yap best alignments...", file=sys.stderr)
    yap_reads = collect_yap_best(args.yap_bam, keys)

    same_strict = 0
    same_chrom_diff_start = 0
    diff_chrom = 0
    missing_bh = 0
    missing_yap = 0
    strand_mismatch_same_start = 0
    same_chrom_same_strand_within_tol = 0

    max_deltas_diff_span: list[int] = []
    start_abs_diff_span: list[int] = []
    end_abs_diff_span: list[int] = []

    detail_rows = []

    for bid, is_r1 in zero_rows:
        key = (bid, is_r1)
        b = bh_loc.get(key)
        y = yap_reads.get(key)
        if b is None:
            missing_bh += 1
            detail_rows.append(
                (
                    bid,
                    int(is_r1),
                    "missing_bhmem",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                )
            )
            continue
        if y is None:
            missing_yap += 1
            bref, bs, be, brv = b
            detail_rows.append(
                (
                    bid,
                    int(is_r1),
                    "missing_yap",
                    bref,
                    bs,
                    be,
                    int(brv),
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                )
            )
            continue
        yloc = (
            y.reference_name,
            y.reference_start,
            y.reference_end,
            y.is_reverse,
        )
        bref, bs, be, brv = b
        yref, ys, ye, yrv = yloc

        if bref == yref:
            d_start = ys - bs
            d_end = ye - be
            max_abs = max(abs(d_start), abs(d_end))
        else:
            d_start = d_end = max_abs = None

        if bref != yref:
            diff_chrom += 1
            status = "diff_chrom"
        elif bs != ys or be != ye:
            same_chrom_diff_start += 1
            status = "same_chrom_diff_span"
            if max_abs is not None:
                max_deltas_diff_span.append(max_abs)
                start_abs_diff_span.append(abs(ys - bs))
                end_abs_diff_span.append(abs(ye - be))
        else:
            if brv != yrv:
                strand_mismatch_same_start += 1
                status = "same_span_diff_strand"
            else:
                same_strict += 1
                status = "same"

        if (
            bref == yref
            and brv == yrv
            and max_abs is not None
            and max_abs <= tol
        ):
            same_chrom_same_strand_within_tol += 1

        ds_out = "" if d_start is None else d_start
        de_out = "" if d_end is None else d_end
        ma_out = "" if max_abs is None else max_abs

        detail_rows.append(
            (
                bid,
                int(is_r1),
                status,
                bref,
                bs,
                be,
                int(brv),
                yref,
                ys,
                ye,
                int(yrv),
                ds_out,
                de_out,
                ma_out,
            )
        )

    n = len(zero_rows)
    print("\n=== Summary (strict: same chrom, start, end, strand) ===")
    print(f"total_zero_zero_reads\t{n}")
    print(f"same_location_strict\t{same_strict}")
    print(f"same_chrom_different_coordinates\t{same_chrom_diff_start}")
    print(f"different_chromosome\t{diff_chrom}")
    print(f"same_span_but_opposite_strand_flag\t{strand_mismatch_same_start}")
    print(f"missing_bhmem_alignment\t{missing_bh}")
    print(f"missing_yap_alignment\t{missing_yap}")

    print(f"\n=== Tolerance (same chrom + same strand, max |Δstart|,|Δend| <= {tol} bp) ===")
    print(f"same_chrom_same_strand_within_{tol}bp\t{same_chrom_same_strand_within_tol}")

    if max_deltas_diff_span:
        print("\n=== Same chrom, different span (status=same_chrom_diff_span) ===")
        print(f"n\t{len(max_deltas_diff_span)}")
        print(
            "# |Δstart| and |Δend| are yap vs bhmem on the same chromosome; "
            "large values = different loci on that chr."
        )
        for label, vals in (
            ("abs_delta_start", start_abs_diff_span),
            ("abs_delta_end", end_abs_diff_span),
            (
                "max_of_abs_start_and_end",
                max_deltas_diff_span,
            ),
        ):
            print(f"{label}_median_bp\t{float(statistics.median(vals)):.2f}")
            print(f"{label}_mean_bp\t{float(statistics.mean(vals)):.2f}")
            print(f"{label}_min_bp\t{min(vals)}")
            print(f"{label}_max_bp\t{max(vals)}")
    else:
        print("\n=== Same chrom, different span: (none) ===")

    if args.output_tsv:
        with open(args.output_tsv, "w", newline="") as fout:
            w = csv.writer(fout, delimiter="\t")
            w.writerow(
                [
                    "base_id",
                    "is_r1",
                    "status",
                    "bhmem_ref",
                    "bhmem_start",
                    "bhmem_end",
                    "bhmem_rev",
                    "yap_ref",
                    "yap_start",
                    "yap_end",
                    "yap_rev",
                    "delta_start_yap_minus_bhmem",
                    "delta_end_yap_minus_bhmem",
                    "max_abs_delta_bp",
                ]
            )
            for row in detail_rows:
                w.writerow(row)
        print(f"Wrote {args.output_tsv}", file=sys.stderr)

    if args.pdf_path:
        write_report_pdf(
            args.pdf_path,
            per_read_tsv=args.per_read_tsv,
            n_total=n,
            same_strict=same_strict,
            same_chrom_diff_span=same_chrom_diff_start,
            diff_chrom=diff_chrom,
            strand_mismatch=strand_mismatch_same_start,
            missing_bh=missing_bh,
            missing_yap=missing_yap,
            tol=tol,
            within_tol=same_chrom_same_strand_within_tol,
            max_deltas_diff_span=max_deltas_diff_span,
        )
        print(f"Wrote {args.pdf_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
