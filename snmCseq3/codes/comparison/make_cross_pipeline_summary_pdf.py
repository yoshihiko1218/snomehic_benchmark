#!/usr/bin/env python3
"""Generate summary PDF for bhmem vs yap cross-pipeline comparison.

Bisulfite-aware NM for both pipelines:
- bhmem: recomputed from BAM read vs unconverted reference, skipping
  C(ref)->T(read) on forward strand, G(ref)->A(read) on reverse strand.
- yap: NM:i tag minus lowercase bisulfite positions in XM tag (h, x, z).

Both methods give identical results for reads at the same position+CIGAR (99.7%).
"""

import os
import sys

import numpy as np
import pysam

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3"
GENOME = "/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa"


def normalize_yap(qname):
    parts = qname.split("_")
    if len(parts) < 2:
        return None
    return (parts[0], parts[1] == "2")


def yap_bs_nm(r):
    """Bisulfite-aware NM from yap/Bismark BAM: NM - lowercase XM positions."""
    nm = int(r.get_tag("NM")) if r.has_tag("NM") else -1
    if nm < 0:
        return -1
    xm = str(r.get_tag("XM")) if r.has_tag("XM") else ""
    converted = sum(1 for c in xm if c in "hxz")
    return nm - converted


def bhmem_bs_nm(r, genome_fa):
    """Bisulfite-aware NM from bhmem BAM: recompute against unconverted ref,
    skipping C>T (forward) and G>A (reverse) bisulfite mismatches."""
    seq = r.query_sequence.upper()
    pairs = r.get_aligned_pairs(matches_only=True)
    ref_seq = genome_fa.fetch(
        r.reference_name, r.reference_start, r.reference_end
    ).upper()
    nm = 0
    for qpos, rpos in pairs:
        rb = ref_seq[rpos - r.reference_start]
        qb = seq[qpos]
        if qb == rb:
            continue
        if not r.is_reverse and rb == "C" and qb == "T":
            continue
        if r.is_reverse and rb == "G" and qb == "A":
            continue
        nm += 1
    nm += sum(l for op, l in r.cigartuples if op in (1, 2))
    return nm


def load_data(bhmem_bam, yap_bam, genome_path, max_reads=20000):
    genome_fa = pysam.FastaFile(genome_path)

    # Scan yap
    print("Scanning yap...", flush=True)
    bam_y = pysam.AlignmentFile(yap_bam, "rb")
    yap = {}
    for r in bam_y:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        parsed = normalize_yap(r.query_name)
        if parsed is None:
            continue
        nm_bs = yap_bs_nm(r)
        if nm_bs < 0:
            continue
        entry = {
            "nm": nm_bs,
            "chrom": r.reference_name,
            "pos": r.reference_start,
            "cigar": r.cigarstring,
            "mapq": int(r.mapping_quality),
        }
        if parsed not in yap or nm_bs < yap[parsed]["nm"]:
            yap[parsed] = entry
    bam_y.close()
    print("  %d yap reads" % len(yap), flush=True)

    # Scan bhmem
    print("Scanning bhmem...", flush=True)
    bam_b = pysam.AlignmentFile(bhmem_bam, "rb")
    results = []
    n = 0
    for r in bam_b:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if not r.has_tag("NM"):
            continue
        is_r2 = bool(r.is_paired and r.is_read2)
        key = (r.query_name, is_r2)
        if key not in yap:
            continue
        yv = yap[key]
        same_pos = r.reference_name == yv["chrom"] and r.reference_start == yv["pos"]
        same_cigar = same_pos and r.cigarstring == yv["cigar"]
        bnm = bhmem_bs_nm(r, genome_fa)
        results.append({
            "is_r2": is_r2,
            "bhmem_nm": bnm,
            "yap_nm": yv["nm"],
            "bhmem_mapq": int(r.mapping_quality),
            "yap_mapq": yv["mapq"],
            "same_pos": same_pos,
            "same_cigar": same_cigar,
        })
        n += 1
        if n % 10000 == 0:
            print("  %d reads..." % n, flush=True)
        if n >= max_reads:
            break
    bam_b.close()
    genome_fa.close()
    print("  %d shared reads" % len(results), flush=True)
    return results


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def make_mapq_2x2_table(ax, data, title):
    ax.axis("off")
    n = len(data)
    if n == 0:
        ax.set_title(title + ": no data")
        return
    hh = sum(1 for d in data if d["bhmem_mapq"] >= 30 and d["yap_mapq"] >= 30)
    hl = sum(1 for d in data if d["bhmem_mapq"] >= 30 and d["yap_mapq"] < 30)
    lh = sum(1 for d in data if d["bhmem_mapq"] < 30 and d["yap_mapq"] >= 30)
    ll = sum(1 for d in data if d["bhmem_mapq"] < 30 and d["yap_mapq"] < 30)
    cell_text = [
        ["%d (%.1f%%)" % (hh, 100.0 * hh / n), "%d (%.1f%%)" % (hl, 100.0 * hl / n)],
        ["%d (%.1f%%)" % (lh, 100.0 * lh / n), "%d (%.1f%%)" % (ll, 100.0 * ll / n)],
    ]
    colors = [["#c6efce", "#ffc7ce"], ["#ffc7ce", "#ffffcc"]]
    table = ax.table(
        cellText=cell_text,
        rowLabels=["bhmem MQ>=30", "bhmem MQ<30"],
        colLabels=["yap MQ>=30", "yap MQ<30"],
        cellColours=colors, loc="center", cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)
    ax.set_title("%s (n=%d)" % (title, n), fontsize=13, fontweight="bold", pad=20)


def make_location_bar(ax, data, title):
    n = len(data)
    if n == 0:
        ax.set_title(title + ": no data")
        return
    sc = sum(1 for d in data if d["same_cigar"])
    sp = sum(1 for d in data if d["same_pos"]) - sc
    dp = n - sc - sp
    colors = ["#4472c4", "#ed7d31", "#a5a5a5"]
    ax.barh([0], [sc], color=colors[0])
    ax.barh([0], [sp], left=[sc], color=colors[1])
    ax.barh([0], [dp], left=[sc + sp], color=colors[2])
    ax.set_xlim(0, n)
    ax.set_yticks([])
    ax.set_xlabel("Number of reads")
    ax.set_title("%s (n=%d)" % (title, n), fontsize=13, fontweight="bold")
    bars = [sc, sp, dp]
    labels = [
        "Same CIGAR\n%.1f%%" % (100.0 * sc / n),
        "Same pos\ndiff CIGAR\n%.1f%%" % (100.0 * sp / n),
        "Diff pos\n%.1f%%" % (100.0 * dp / n),
    ]
    positions = [sc / 2, sc + sp / 2, sc + sp + dp / 2]
    for pos, lbl, count in zip(positions, labels, bars):
        if count > n * 0.05:
            ax.text(pos, 0, lbl, ha="center", va="center", fontsize=9, fontweight="bold")


def make_nm_scatter(ax, data, title, box_color):
    if not data:
        ax.set_title(title + ": no data")
        return
    b = np.array([d["bhmem_nm"] for d in data], dtype=float)
    y = np.array([d["yap_nm"] for d in data], dtype=float)
    n = len(data)
    bl = int(np.sum(b < y))
    eq = int(np.sum(b == y))
    yl = int(np.sum(y < b))
    valid = n > 1 and np.std(b) > 0 and np.std(y) > 0
    pearson = float(np.corrcoef(b, y)[0, 1]) if valid else float("nan")
    ax.scatter(y, b, alpha=0.12, s=6, rasterized=True, edgecolors="none")
    mx = max(b.max(), y.max(), 1)
    ax.plot([0, mx], [0, mx], "k--", alpha=0.4, linewidth=1)
    ax.set_xlabel("yap BS-aware NM")
    ax.set_ylabel("bhmem BS-aware NM")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlim(-0.5, min(mx, 50) + 0.5)
    ax.set_ylim(-0.5, min(mx, 50) + 0.5)
    stats = (
        "n = %d\n"
        "bhmem < yap: %.1f%%\n"
        "equal: %.1f%%\n"
        "yap < bhmem: %.1f%%\n"
        "mean bhmem: %.2f\n"
        "mean yap: %.2f\n"
        "Pearson r: %.4f"
    ) % (n, 100.0 * bl / n, 100.0 * eq / n, 100.0 * yl / n,
         np.mean(b), np.mean(y), pearson)
    ax.text(
        0.02, -0.20, stats, transform=ax.transAxes,
        fontsize=9, verticalalignment="top", fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor=box_color, alpha=0.8),
    )


def make_mapq_scatter(ax, data, title):
    if not data:
        ax.set_title(title + ": no data")
        return
    b = np.array([d["bhmem_mapq"] for d in data], dtype=float)
    y = np.array([d["yap_mapq"] for d in data], dtype=float)
    n = len(data)
    jitter = 0.5
    bj = b + np.random.uniform(-jitter, jitter, n)
    yj = y + np.random.uniform(-jitter, jitter, n)
    ax.scatter(yj, bj, alpha=0.08, s=4, rasterized=True, edgecolors="none")
    ax.plot([0, 60], [0, 60], "k--", alpha=0.4, linewidth=1)
    ax.axhline(30, color="red", alpha=0.3, linewidth=0.8, linestyle=":")
    ax.axvline(30, color="red", alpha=0.3, linewidth=0.8, linestyle=":")
    ax.set_xlabel("yap MAPQ")
    ax.set_ylabel("bhmem MAPQ")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlim(-1, 62)
    ax.set_ylim(-1, 62)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bhmem-bam",
                    default=os.path.join(BASE, "04.bhmem_bam/SRR21549292.bhmem.bam"))
    ap.add_argument("--yap-bam",
                    default=os.path.join(BASE, "alignment/Group22/bam/SRR21549292.3C.sorted.bam"))
    ap.add_argument("--genome", default=GENOME)
    ap.add_argument("-o", "--output",
                    default=os.path.join(BASE, "mapq_comparison/SRR21549292/cross_pipeline_summary.pdf"))
    ap.add_argument("--max-reads", type=int, default=20000)
    args = ap.parse_args()

    data = load_data(args.bhmem_bam, args.yap_bam, args.genome, args.max_reads)
    r1 = [d for d in data if not d["is_r2"]]
    r2 = [d for d in data if d["is_r2"]]

    print("Generating PDF (%d reads)..." % len(data), flush=True)

    with PdfPages(args.output) as pdf:
        # Page 1: Location concordance + MAPQ 2x2
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        for col, (label, sub) in enumerate([("R1", r1), ("R2", r2), ("All", data)]):
            make_location_bar(axes[0, col], sub, label)
        for col, (label, sub) in enumerate([("R1", r1), ("R2", r2), ("All", data)]):
            make_mapq_2x2_table(axes[1, col], sub, label)
        fig.suptitle("Location concordance & MAPQ agreement\n"
                     "bhmem vs yap (shared reads, n=%d)" % len(data),
                     fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 2: BS-aware NM scatter — all shared reads
        fig, axes = plt.subplots(1, 3, figsize=(18, 8))
        for col, (label, sub) in enumerate([("R1", r1), ("R2", r2), ("All", data)]):
            make_nm_scatter(axes[col], sub, label, "lightcyan")
        fig.suptitle("Bisulfite-aware NM — all shared reads\n"
                     "bhmem: recomputed vs unconverted ref | yap: NM - lowercase XM",
                     fontsize=13, fontweight="bold")
        fig.subplots_adjust(bottom=0.22)
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 3: BS-aware NM scatter — same CIGAR only
        sc_r1 = [d for d in r1 if d["same_cigar"]]
        sc_r2 = [d for d in r2 if d["same_cigar"]]
        sc_all = sc_r1 + sc_r2
        fig, axes = plt.subplots(1, 3, figsize=(18, 8))
        for col, (label, sub) in enumerate([("R1", sc_r1), ("R2", sc_r2), ("All", sc_all)]):
            make_nm_scatter(axes[col], sub, label, "wheat")
        fig.suptitle("Bisulfite-aware NM — same position + same CIGAR only\n"
                     "Reads placed identically by both pipelines",
                     fontsize=13, fontweight="bold")
        fig.subplots_adjust(bottom=0.22)
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 4: MAPQ scatter
        fig, axes = plt.subplots(1, 3, figsize=(18, 7))
        for col, (label, sub) in enumerate([("R1", r1), ("R2", r2), ("All", data)]):
            make_mapq_scatter(axes[col], sub, label)
        fig.suptitle("MAPQ comparison — all shared reads\n"
                     "Red dotted lines at MAPQ=30",
                     fontsize=13, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 5: BS-aware NM — discrepant MAPQ (bhmem<30, yap>=30)
        disc_r1 = [d for d in r1 if d["bhmem_mapq"] < 30 and d["yap_mapq"] >= 30]
        disc_r2 = [d for d in r2 if d["bhmem_mapq"] < 30 and d["yap_mapq"] >= 30]
        disc_all = disc_r1 + disc_r2
        fig, axes = plt.subplots(1, 3, figsize=(18, 8))
        for col, (label, sub) in enumerate([("R1", disc_r1), ("R2", disc_r2), ("All", disc_all)]):
            make_nm_scatter(axes[col], sub, label, "mistyrose")
        fig.suptitle("Bisulfite-aware NM — discrepant MAPQ (bhmem MQ<30, yap MQ>=30)\n"
                     "Reads where yap is confident but bhmem is not",
                     fontsize=13, fontweight="bold")
        fig.subplots_adjust(bottom=0.22)
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 6: Summary tables
        fig = plt.figure(figsize=(18, 22))

        # --- Location concordance table ---
        ax1 = fig.add_axes([0.05, 0.78, 0.9, 0.12])
        ax1.axis("off")
        ax1.set_title("Table 1: Location concordance", fontsize=14, fontweight="bold",
                       loc="left", pad=10)
        loc_rows = []
        for label, sub in [("R1", r1), ("R2", r2), ("All", data)]:
            n = len(sub)
            sc = sum(1 for d in sub if d["same_cigar"])
            sp = sum(1 for d in sub if d["same_pos"]) - sc
            dp = n - sc - sp
            loc_rows.append([
                label, str(n),
                "%d (%.1f%%)" % (sc, 100.0 * sc / n),
                "%d (%.1f%%)" % (sp, 100.0 * sp / n),
                "%d (%.1f%%)" % (dp, 100.0 * dp / n),
            ])
        t1 = ax1.table(
            cellText=loc_rows,
            colLabels=["Mate", "Total", "Same pos + CIGAR", "Same pos, diff CIGAR", "Diff position"],
            loc="center", cellLoc="center",
        )
        t1.auto_set_font_size(False)
        t1.set_fontsize(10)
        t1.scale(1, 1.6)

        # --- NM comparison table (all reads) ---
        ax2 = fig.add_axes([0.05, 0.62, 0.9, 0.12])
        ax2.axis("off")
        ax2.set_title("Table 2: BS-aware NM comparison (all shared reads)", fontsize=14,
                       fontweight="bold", loc="left", pad=10)
        nm_rows = []
        for label, sub in [("R1", r1), ("R2", r2), ("All", data)]:
            n = len(sub)
            bl = sum(1 for d in sub if d["bhmem_nm"] < d["yap_nm"])
            eq = sum(1 for d in sub if d["bhmem_nm"] == d["yap_nm"])
            yl = sum(1 for d in sub if d["yap_nm"] < d["bhmem_nm"])
            nm_rows.append([
                label, str(n),
                "%d (%.1f%%)" % (bl, 100.0 * bl / n),
                "%d (%.1f%%)" % (eq, 100.0 * eq / n),
                "%d (%.1f%%)" % (yl, 100.0 * yl / n),
            ])
        t2 = ax2.table(
            cellText=nm_rows,
            colLabels=["Mate", "Total", "bhmem < yap", "Equal", "yap < bhmem"],
            loc="center", cellLoc="center",
        )
        t2.auto_set_font_size(False)
        t2.set_fontsize(10)
        t2.scale(1, 1.6)

        # --- NM comparison table (same CIGAR only) ---
        ax3 = fig.add_axes([0.05, 0.46, 0.9, 0.12])
        ax3.axis("off")
        ax3.set_title("Table 3: BS-aware NM comparison (same pos + same CIGAR only)",
                       fontsize=14, fontweight="bold", loc="left", pad=10)
        nm_sc_rows = []
        for label, sub in [("R1", sc_r1), ("R2", sc_r2), ("All", sc_all)]:
            n = len(sub)
            if n == 0:
                nm_sc_rows.append([label, "0", "-", "-", "-"])
                continue
            bl = sum(1 for d in sub if d["bhmem_nm"] < d["yap_nm"])
            eq = sum(1 for d in sub if d["bhmem_nm"] == d["yap_nm"])
            yl = sum(1 for d in sub if d["yap_nm"] < d["bhmem_nm"])
            nm_sc_rows.append([
                label, str(n),
                "%d (%.1f%%)" % (bl, 100.0 * bl / n),
                "%d (%.1f%%)" % (eq, 100.0 * eq / n),
                "%d (%.1f%%)" % (yl, 100.0 * yl / n),
            ])
        t3 = ax3.table(
            cellText=nm_sc_rows,
            colLabels=["Mate", "Total", "bhmem < yap", "Equal", "yap < bhmem"],
            loc="center", cellLoc="center",
        )
        t3.auto_set_font_size(False)
        t3.set_fontsize(10)
        t3.scale(1, 1.6)

        # --- MAPQ 2x2 tables ---
        ax4 = fig.add_axes([0.05, 0.28, 0.9, 0.14])
        ax4.axis("off")
        ax4.set_title("Table 4: MAPQ 2x2 (MAPQ >= 30 vs < 30)", fontsize=14,
                       fontweight="bold", loc="left", pad=10)
        mapq_rows = []
        for label, sub in [("R1", r1), ("R2", r2), ("All", data)]:
            n = len(sub)
            hh = sum(1 for d in sub if d["bhmem_mapq"] >= 30 and d["yap_mapq"] >= 30)
            hl = sum(1 for d in sub if d["bhmem_mapq"] >= 30 and d["yap_mapq"] < 30)
            lh = sum(1 for d in sub if d["bhmem_mapq"] < 30 and d["yap_mapq"] >= 30)
            ll = sum(1 for d in sub if d["bhmem_mapq"] < 30 and d["yap_mapq"] < 30)
            mapq_rows.append([
                label, str(n),
                "%d (%.1f%%)" % (hh, 100.0 * hh / n),
                "%d (%.1f%%)" % (hl, 100.0 * hl / n),
                "%d (%.1f%%)" % (lh, 100.0 * lh / n),
                "%d (%.1f%%)" % (ll, 100.0 * ll / n),
            ])
        t4 = ax4.table(
            cellText=mapq_rows,
            colLabels=["Mate", "Total", "Both >= 30", "bhmem>=30\nyap<30",
                        "bhmem<30\nyap>=30", "Both < 30"],
            loc="center", cellLoc="center",
        )
        t4.auto_set_font_size(False)
        t4.set_fontsize(10)
        t4.scale(1, 1.6)

        fig.suptitle("Summary tables — bhmem vs yap (n=%d shared reads)" % len(data),
                     fontsize=16, fontweight="bold", y=0.95)
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 7: Methods — NM correction explanation
        fig = plt.figure(figsize=(18, 14))
        ax = fig.add_axes([0.08, 0.05, 0.84, 0.85])
        ax.axis("off")

        methods_text = """Bisulfite-aware NM correction method

Problem
  Both bhmem (BWA-meth) and yap (Bismark/Bowtie2) align bisulfite-converted queries
  against converted reference genomes. The NM:i tag in each BAM has a different meaning:

    - Bismark NM: mismatches against the unconverted genome (includes all bisulfite C>T)
    - Bhmem NM:   mismatches against the converted genome used by BWA
                   (excludes one direction of bisulfite but includes the other)

  Neither NM tag is directly comparable across pipelines.

Solution: Unified bisulfite-aware NM
  For each read, recompute NM against the unconverted reference genome (mm10.fa),
  skipping mismatches that are expected bisulfite conversions:

    - Forward strand (is_reverse = False): skip C(ref) > T(read)
    - Reverse strand (is_reverse = True):  skip G(ref) > A(read)

  This rule applies identically to both R1 and R2, using the BAM strand flag as-is.

  For yap/Bismark reads, an equivalent shortcut exists:
    BS-aware NM = NM:i tag  -  count(lowercase h, x, z in XM tag)
  where lowercase XM characters mark unmethylated cytosine positions (bisulfite-converted).

Validation
  On reads placed at the same position with the same CIGAR by both pipelines:
    - R1: 99.7% exact NM agreement between the two methods
    - R2: 100% exact NM agreement
  The two corrections are mathematically equivalent when the alignment is identical.

Implementation
  bhmem_bs_nm(read, genome_fa):
    For each aligned (query_pos, ref_pos) pair from CIGAR:
      ref_base = genome_fa[chrom][ref_pos]
      read_base = read.query_sequence[query_pos]
      if ref_base == read_base: continue          # match
      if forward and ref=C and read=T: continue   # bisulfite
      if reverse and ref=G and read=A: continue   # bisulfite (complement)
      nm += 1                                     # real mismatch
    nm += indel_bases_from_CIGAR
    return nm

  yap_bs_nm(read):
    return read.NM - count(c for c in read.XM if c in "hxz")

Note on bhmem PBAT strand handling
  Bhmem with -pbat flips the strand flag for R1 reads (Bhmem.java line 650).
  However, empirical testing confirms no additional flip is needed for the NM
  correction formula above -- the BAM strand flag can be used directly for
  both R1 and R2.
"""
        ax.text(0.02, 0.98, methods_text, transform=ax.transAxes,
                fontsize=11, verticalalignment="top", fontfamily="monospace",
                linespacing=1.4)
        fig.suptitle("Methods: Bisulfite-aware NM correction",
                     fontsize=16, fontweight="bold")
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

    print("Wrote %s" % args.output)


if __name__ == "__main__":
    main()
