#!/usr/bin/env python3
"""
For reads with bhmem MAPQ < 30 and yap MAPQ > 30 (one cell):
  - NM from bhmem (aligner tag)
  - Bisulfite-corrected substitution mismatch on bhmem (CIGAR M/=/X only; see bisulfite_corrected_mismatch.py)
  - raw + bisulfite-corrected mismatch counts on yap (aligned_pairs / matches_only)

Outputs: per-read TSV, summary TXT, PDF with histograms and scatter/hexbin.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import numpy as np
    import pysam
except ImportError as e:
    print("ERROR: need pysam and numpy:", e, file=sys.stderr)
    sys.exit(1)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    plt = None
    PdfPages = None

from bisulfite_corrected_mismatch import count_mismatches_corrected_cigar_mx


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


def count_mismatches_corrected(read, fasta):
    """
    YAP/Bismark bisulfite-corrected mismatch heuristic.

    Prefer Bismark's conversion class tag `XR:Z`:
      - `XR:Z:CT` -> ignore ref C vs read T substitutions
      - `XR:Z:GA` -> ignore ref G vs read A substitutions

    If `XR` is missing, fall back to the older mapping-strand heuristic.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0, 0, 0

    chrom = read.reference_name
    raw_mm = corr_mm = bs_ign = 0
    rev = read.is_reverse

    xr = ""
    if read.has_tag("XR"):
        xr = str(read.get_tag("XR")).upper()

    # XR typically is one of {'CT','GA'} for bismark bisulfite conversion class.
    ignore_ct = xr == "CT"
    ignore_ga = xr == "GA"

    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if qpos is None or rpos is None:
            continue
        rb = read.query_sequence[qpos].upper()
        refb = fasta.fetch(chrom, rpos, rpos + 1).upper()
        if rb == refb:
            continue

        raw_mm += 1

        is_conv = False
        if ignore_ct:
            is_conv = refb == "C" and rb == "T"
        elif ignore_ga:
            is_conv = refb == "G" and rb == "A"
        else:
            # Fallback: original strand-based rule.
            if not rev:
                is_conv = refb == "C" and rb == "T"
            else:
                is_conv = refb == "G" and rb == "A"

        if is_conv:
            bs_ign += 1
        else:
            corr_mm += 1

    return raw_mm, corr_mm, bs_ign


def get_nm(read) -> int | None:
    if not read.has_tag("NM"):
        return None
    return int(read.get_tag("NM"))


def load_subset(path: str):
    rows = []
    keys = set()
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            bid = row["base_id"]
            is_r1 = bool(int(row["is_r1"]))
            mb = int(row["mapq_bhmem"])
            my = int(row["mapq_yap"])
            keys.add((bid, is_r1))
            rows.append(
                {
                    "base_id": bid,
                    "is_r1": is_r1,
                    "mapq_bhmem": mb,
                    "mapq_yap": my,
                }
            )
    return rows, keys


def collect_bhmem_primary(bam_path: str, keys: set) -> dict:
    """key -> primary AlignedSegment (same selection as former NM collection)."""
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            is_r1 = bool(read.flag & 64)
            if not (read.flag & 64) and not (read.flag & 128):
                continue
            key = (qname, is_r1)
            if key not in keys:
                continue
            out[key] = read
    return out


def collect_yap_best(bam_path: str, keys: set) -> dict:
    """key -> best AlignedSegment (max MAPQ, prefer primary no -l/-r/-m)."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "yap_high_bhmem_low_tsv",
        help="yap_high_bhmem_low.tsv (base_id, is_r1, mapq_bhmem, mapq_yap, ...)",
    )
    ap.add_argument("bhmem_bam", help="bhmem BAM")
    ap.add_argument("yap_bam", help="yap 3C sorted BAM (Bowtie2)")
    ap.add_argument("reference_fasta", help="Genomic FASTA (indexed)")
    ap.add_argument("-o", "--output-prefix", default="discrepant_mismatch_report")
    ap.add_argument("--max-reads", type=int, default=0, help="Limit rows for testing (0=all)")
    args = ap.parse_args()

    rows, keys = load_subset(args.yap_high_bhmem_low_tsv)
    if args.max_reads:
        rows = rows[: args.max_reads]
        keys = {(r["base_id"], r["is_r1"]) for r in rows}

    progress_every = 100000
    if args.max_reads and args.max_reads <= 200000:
        progress_every = 1000

    print(f"Loaded {len(rows)} discrepant reads", file=sys.stderr)

    print("Collecting bhmem primary alignments...", file=sys.stderr)
    bhmem_read = collect_bhmem_primary(args.bhmem_bam, keys)
    print(f"  bhmem keys found: {len(bhmem_read)}", file=sys.stderr)

    print("Collecting yap alignments (best per read)...", file=sys.stderr)
    yap_best = collect_yap_best(args.yap_bam, keys)
    print(f"  yap keys found: {len(yap_best)}", file=sys.stderr)

    fasta = pysam.FastaFile(args.reference_fasta)

    per_read_path = f"{args.output_prefix}.per_read.tsv"
    nm_bh = []
    raw_bh = []
    corr_bh = []
    raw_y = []
    corr_y = []
    mqb = []
    mqy = []

    with open(per_read_path, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(
            [
                "base_id",
                "is_r1",
                "mapq_bhmem",
                "mapq_yap",
                "nm_bhmem",
                "bhmem_raw_mismatch",
                "bhmem_corrected_mismatch",
                "bhmem_bisulfite_ignored",
                "yap_raw_mismatch",
                "yap_corrected_mismatch",
                "yap_bisulfite_ignored",
            ]
        )
        n_done = 0
        for rec in rows:
            key = (rec["base_id"], rec["is_r1"])
            bread = bhmem_read.get(key)
            yread = yap_best.get(key)
            if bread is None or yread is None:
                w.writerow(
                    [
                        rec["base_id"],
                        int(rec["is_r1"]),
                        rec["mapq_bhmem"],
                        rec["mapq_yap"],
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                    ]
                )
                continue
            nmb = get_nm(bread)
            if nmb is None:
                nmb = ""
            braw, bcorr, bbs = count_mismatches_corrected_cigar_mx(bread, fasta)
            raw, corr, bs = count_mismatches_corrected(yread, fasta)
            w.writerow(
                [
                    rec["base_id"],
                    int(rec["is_r1"]),
                    rec["mapq_bhmem"],
                    rec["mapq_yap"],
                    nmb,
                    braw,
                    bcorr,
                    bbs,
                    raw,
                    corr,
                    bs,
                ]
            )
            nm_bh.append(float(nmb) if nmb != "" else np.nan)
            raw_bh.append(braw)
            corr_bh.append(bcorr)
            raw_y.append(raw)
            corr_y.append(corr)
            mqb.append(rec["mapq_bhmem"])
            mqy.append(rec["mapq_yap"])
            n_done += 1
            if progress_every and n_done % progress_every == 0:
                print(f"  processed {n_done} / {len(rows)}...", file=sys.stderr)

    fasta.close()

    nm_bh = np.array(nm_bh, dtype=np.float64)
    raw_bh = np.array(raw_bh, dtype=np.float64)
    corr_bh = np.array(corr_bh, dtype=np.float64)
    raw_y = np.array(raw_y, dtype=np.float64)
    corr_y = np.array(corr_y, dtype=np.float64)

    def _corr(a: np.ndarray, b: np.ndarray) -> float:
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 2:
            return float("nan")
        return float(np.corrcoef(a[m], b[m])[0, 1])

    summary_path = f"{args.output_prefix}.summary.txt"
    with open(summary_path, "w") as sf:
        sf.write(f"n_reads_with_both_alignments\t{len(nm_bh)}\n")
        sf.write(f"nm_bhmem_mean\t{float(np.nanmean(nm_bh)):.6f}\n")
        sf.write(f"nm_bhmem_median\t{float(np.nanmedian(nm_bh)):.6f}\n")
        sf.write(f"nm_bhmem_std\t{float(np.nanstd(nm_bh)):.6f}\n")
        sf.write(f"bhmem_corrected_mismatch_mean\t{float(np.mean(corr_bh)):.6f}\n")
        sf.write(f"bhmem_corrected_mismatch_median\t{float(np.median(corr_bh)):.6f}\n")
        sf.write(f"bhmem_corrected_mismatch_std\t{float(np.std(corr_bh)):.6f}\n")
        sf.write(f"yap_raw_mismatch_mean\t{float(np.mean(raw_y)):.6f}\n")
        sf.write(f"yap_raw_mismatch_median\t{float(np.median(raw_y)):.6f}\n")
        sf.write(f"yap_corrected_mismatch_mean\t{float(np.mean(corr_y)):.6f}\n")
        sf.write(f"yap_corrected_mismatch_median\t{float(np.median(corr_y)):.6f}\n")
        sf.write(f"yap_corrected_mismatch_std\t{float(np.std(corr_y)):.6f}\n")
        sf.write(
            f"correlation_nm_bhmem_vs_yap_corrected\t{_corr(nm_bh, corr_y):.6f}\n"
        )
        sf.write(f"correlation_nm_bhmem_vs_yap_raw\t{_corr(nm_bh, raw_y):.6f}\n")
        sf.write(
            f"correlation_bhmem_corrected_vs_yap_corrected\t{_corr(corr_bh, corr_y):.6f}\n"
        )
        sf.write(
            f"correlation_nm_bhmem_vs_bhmem_raw_mismatch\t{_corr(nm_bh, raw_bh):.6f}\n"
        )
        sf.write(
            f"correlation_nm_bhmem_vs_bhmem_corrected_mismatch\t{_corr(nm_bh, corr_bh):.6f}\n"
        )

    print(f"Wrote {per_read_path}", file=sys.stderr)
    print(f"Wrote {summary_path}", file=sys.stderr)

    if plt is None or PdfPages is None:
        print("WARNING: matplotlib missing, skip plots", file=sys.stderr)
        return

    pdf_path = f"{args.output_prefix}.plots.pdf"
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))

    ax = axes[0, 0]
    nm_ok = nm_bh[np.isfinite(nm_bh)]
    if len(nm_ok):
        ax.hist(
            nm_ok,
            bins=min(50, int(np.nanmax(nm_ok)) + 2),
            color="steelblue",
            alpha=0.8,
            edgecolor="black",
        )
    ax.set_xlabel("NM (bhmem)")
    ax.set_ylabel("Count")
    ax.set_title("Bhmem NM tag")

    ax = axes[0, 1]
    ax.hist(
        corr_bh,
        bins=min(40, int(corr_bh.max()) + 2),
        color="#6d4c41",
        alpha=0.8,
        edgecolor="black",
    )
    ax.set_xlabel("Bisulfite-corrected mismatch (bhmem, CIGAR M/=/X)")
    ax.set_ylabel("Count")
    ax.set_title("Bhmem corrected mismatch")

    ax = axes[0, 2]
    ax.hist(raw_y, bins=min(50, int(raw_y.max()) + 2), color="coral", alpha=0.8, edgecolor="black")
    ax.set_xlabel("Raw mismatch vs genome (yap)")
    ax.set_ylabel("Count")
    ax.set_title("Yap raw mismatch (genomic)")

    ax = axes[1, 0]
    ax.hist(corr_y, bins=min(40, int(corr_y.max()) + 2), color="seagreen", alpha=0.8, edgecolor="black")
    ax.set_xlabel("Bisulfite-corrected mismatch (yap)")
    ax.set_ylabel("Count")
    ax.set_title("Yap corrected mismatch")

    ax = axes[1, 1]
    m2 = np.isfinite(nm_bh)
    hb = ax.hexbin(nm_bh[m2], corr_y[m2], gridsize=40, cmap="viridis", mincnt=1)
    ax.set_xlabel("NM bhmem")
    ax.set_ylabel("Yap corrected mismatch")
    ax.set_title("NM (bhmem) vs yap corrected mismatch")
    plt.colorbar(hb, ax=ax, label="Count")

    ax = axes[1, 2]
    hb2 = ax.hexbin(corr_bh, corr_y, gridsize=40, cmap="plasma", mincnt=1)
    ax.set_xlabel("Bhmem corrected mismatch")
    ax.set_ylabel("Yap corrected mismatch")
    ax.set_title("Bhmem vs yap bisulfite-corrected (substitutions)")
    plt.colorbar(hb2, ax=ax, label="Count")

    plt.suptitle(
        "Discrepant reads: bhmem MAPQ<30, yap MAPQ>30",
        fontsize=11,
    )
    plt.tight_layout()

    pr_nm_corr = _corr(nm_bh, corr_bh)
    fig2, ax_s = plt.subplots(figsize=(6.5, 6))
    mplot = np.isfinite(nm_bh)
    ax_s.scatter(
        nm_bh[mplot],
        corr_bh[mplot],
        alpha=0.22,
        s=12,
        c="#1a365d",
        edgecolors="none",
        rasterized=True,
    )
    ax_s.set_xlabel("NM (bhmem)")
    ax_s.set_ylabel("Bhmem bisulfite-corrected mismatch (CIGAR M/=/X)")
    ax_s.set_title(
        "Pearson correlation: NM vs bhmem corrected mismatch\n"
        f"r = {pr_nm_corr:.4f}  (n = {int(mplot.sum())})"
    )
    mx = float(max(np.nanmax(nm_bh), float(np.max(corr_bh))))
    ax_s.plot([0, mx], [0, mx], "k--", alpha=0.35, lw=1, label="y = x")
    ax_s.set_xlim(left=-0.5)
    ax_s.set_ylim(bottom=-0.5)
    ax_s.legend(loc="upper left", fontsize=8)
    fig2.tight_layout()

    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        pdf.savefig(fig2, dpi=150)
        plt.close(fig2)
    print(f"Wrote {pdf_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
