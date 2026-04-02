#!/usr/bin/env python3
"""
Compare **min(converted-genome trials)** NM metric against Bhmem NM tag,
then cross-pipeline (bhmem vs yap) on shared reads.

Section 1: min(trials) vs Bhmem NM:i — agreement rates, scatter plot.
Section 2: min(trials) on bhmem BAM vs min(trials) on yap BAM for the same reads.

Requires bisulfite-converted FASTAs (CT/GA from Bismark Bisulfite_Genome).

Example::

  python report_min_trial_nm.py \\
    /path/to/bhmem.bam \\
    /path/to/yap.3C.sorted.bam \\
    /path/to/Bisulfite_Genome \\
    --joined-tsv mapq_comparison/.../mapq_comparison.joined.tsv \\
    -o output_prefix \\
    --max-reads 20000
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam
from bisulfite_corrected_mismatch import (
    pbat_converted_genome_trial_distances,
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
    _pbat_converted_query_variants,
    _reverse_complement_dna,
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

DEFAULT_BISULFITE = (
    "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"
)


def _normalize_yap_key(qname: str):
    """Normalize yap QNAME to (base_id, is_r2).

    Yap names: ``SRR21549292.1_1_1`` (R1) or ``SRR21549292.1_2_2`` (R2).
    Split reads have suffix like ``-l``, ``-r``, ``-m``.
    Returns ``(base_id, is_r2)`` or ``None`` if unparseable.
    """
    parts = qname.split("_")
    if len(parts) < 2:
        return None
    base = parts[0]
    mate = parts[1]
    if mate == "1":
        return (base, False)
    if mate == "2":
        return (base, True)
    return None


def _all_trial_distances(read, fasta_ct, fasta_ga):
    """Compute edit distances for all conversion × genome trials.

    Unlike ``pbat_converted_genome_trial_distances`` which uses ``read.is_read2``
    to determine query conversion, this tries **both** G→A and C→T conversions
    × both CT and GA genomes × both orientations = up to 8 trials.
    Returns list of valid distances. Suitable for yap SE reads where is_read2
    is always False.
    """
    if read.is_unmapped or read.query_sequence is None:
        return []

    qs = read.query_sequence.upper()
    c_ct = bisulfite_converted_contig_name(fasta_ct, read.reference_name, "CT")
    c_ga = bisulfite_converted_contig_name(fasta_ga, read.reference_name, "GA")
    if c_ct is None or c_ga is None:
        return []

    trials = []
    for conv_fn in (lambda s: s.replace("G", "A"), lambda s: s.replace("C", "T")):
        v1 = conv_fn(qs)
        v2 = _reverse_complement_dna(conv_fn(_reverse_complement_dna(qs)))
        for fa, cname in ((fasta_ct, c_ct), (fasta_ga, c_ga)):
            for qconv in (v1, v2):
                d = count_nm_style_edit_distance_converted_explicit(
                    read, fa, ref_contig=cname, query_converted=qconv
                )
                if d >= 0:
                    trials.append(d)
    return trials


def compute_min_trials(bam_path, ct_fa, ga_fa, read_set=None, max_reads=0,
                       name_mode="bhmem"):
    """Compute min(4 trials) for each primary mapped read.

    Args:
        bam_path: BAM file path.
        ct_fa: pysam.FastaFile for CT-converted genome.
        ga_fa: pysam.FastaFile for GA-converted genome.
        read_set: if not None, only process reads whose (base_id, is_read2)
                  key is in this set. Useful for restricting to shared reads.
        max_reads: stop after this many qualifying reads (0 = no limit).
        name_mode: "bhmem" (QNAME is base_id, flags give R1/R2) or
                   "yap" (QNAME has ``_1``/``_2`` suffix, SE reads).

    Returns:
        dict mapping (base_id, is_read2) -> dict with keys:
            min_trial: int (min of 4 trial distances)
            all_trials: list[int]
            nm_tag: int or None
            mapq: int
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = {}
    n = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None:
            continue

        if name_mode == "yap":
            parsed = _normalize_yap_key(read.query_name)
            if parsed is None:
                continue
            key = parsed
        else:
            is_r2 = bool(read.is_paired and read.is_read2)
            key = (read.query_name, is_r2)

        if read_set is not None and key not in read_set:
            continue

        n += 1
        if max_reads and n > max_reads:
            break
        if n % 5000 == 0:
            print(f"  {bam_path}: {n} reads processed...", flush=True)

        if name_mode == "yap":
            trials = _all_trial_distances(read, ct_fa, ga_fa)
        else:
            trials = pbat_converted_genome_trial_distances(read, ct_fa, ga_fa)
        nm_tag = int(read.get_tag("NM")) if read.has_tag("NM") else None
        mq = int(read.mapping_quality)

        if trials:
            mn = min(trials)
        else:
            mn = nm_tag if nm_tag is not None else -1

        entry = {
            "min_trial": mn,
            "all_trials": trials,
            "nm_tag": nm_tag,
            "mapq": mq,
            "chrom": read.reference_name,
            "pos": read.reference_start,
            "cigar": read.cigarstring,
        }
        # For yap, same base_id can have split reads — keep lowest min_trial
        if key in results:
            if mn < results[key]["min_trial"]:
                results[key] = entry
        else:
            results[key] = entry

    bam.close()
    return results


def load_joined_keys(tsv_path, max_rows=0):
    """Load (base_id, is_r1) pairs from a mapq_comparison joined TSV."""
    keys = set()
    n = 0
    with open(tsv_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            base_id = parts[0]
            is_r1 = int(parts[1])
            is_r2 = not bool(is_r1)
            keys.add((base_id, is_r2))
            n += 1
            if max_rows and n >= max_rows:
                break
    return keys


def section1_bhmem_agreement(results, out_prefix):
    """Section 1: min(trials) vs Bhmem NM:i tag."""
    r1_data = []
    r2_data = []
    for (name, is_r2), v in results.items():
        if v["nm_tag"] is None:
            continue
        entry = (v["min_trial"], v["nm_tag"], v["mapq"])
        if is_r2:
            r2_data.append(entry)
        else:
            r1_data.append(entry)

    lines = []
    lines.append("=" * 60)
    lines.append("Section 1: min(trials) vs Bhmem NM:i tag")
    lines.append("=" * 60)
    lines.append("")

    for label, data in [("R1", r1_data), ("R2", r2_data), ("All", r1_data + r2_data)]:
        if not data:
            continue
        mn_arr = np.array([d[0] for d in data])
        nm_arr = np.array([d[1] for d in data])
        n = len(data)
        exact = int(np.sum(mn_arr == nm_arr))
        mn_le = int(np.sum(mn_arr <= nm_arr))
        mn_gt = int(np.sum(mn_arr > nm_arr))
        diff = mn_arr.astype(float) - nm_arr.astype(float)

        lines.append(f"--- {label} ({n} reads) ---")
        lines.append(f"  Exact match (min==NM):  {exact:7d}  ({100*exact/n:.2f}%)")
        lines.append(f"  min <= NM:              {mn_le:7d}  ({100*mn_le/n:.2f}%)")
        lines.append(f"  min > NM:               {mn_gt:7d}  ({100*mn_gt/n:.2f}%)")
        lines.append(f"  mean(min - NM):         {np.mean(diff):+.3f}")
        lines.append(f"  median(min - NM):       {np.median(diff):+.1f}")
        lines.append(f"  mean |min - NM|:        {np.mean(np.abs(diff)):.3f}")
        if n > 1:
            valid = (np.std(mn_arr) > 0) and (np.std(nm_arr) > 0)
            if valid:
                pearson = float(np.corrcoef(mn_arr, nm_arr)[0, 1])
                lines.append(f"  Pearson r:              {pearson:.4f}")
            else:
                lines.append(f"  Pearson r:              N/A (no variance)")
        lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    summary_path = f"{out_prefix}.section1_summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary_text + "\n")
    print(f"Wrote {summary_path}")

    # TSV
    tsv_path = f"{out_prefix}.section1.per_read.tsv"
    with open(tsv_path, "w") as f:
        f.write("read_name\tis_r2\tmin_trial\tnm_tag\tmapq\tdiff\n")
        for (name, is_r2), v in sorted(results.items()):
            if v["nm_tag"] is None:
                continue
            diff = v["min_trial"] - v["nm_tag"]
            f.write(f"{name}\t{int(is_r2)}\t{v['min_trial']}\t{v['nm_tag']}\t{v['mapq']}\t{diff}\n")
    print(f"Wrote {tsv_path}")

    # Scatter plot
    if HAS_MPL and (r1_data or r2_data):
        fig, axes = plt.subplots(1, 3, figsize=(18, 8))
        for ax, label, data in zip(axes, ["R1", "R2", "All"],
                                    [r1_data, r2_data, r1_data + r2_data]):
            if not data:
                ax.set_title(f"{label}: no data")
                continue
            mn_arr = np.array([d[0] for d in data])
            nm_arr = np.array([d[1] for d in data])
            n = len(data)
            exact = int(np.sum(mn_arr == nm_arr))
            mn_le = int(np.sum(mn_arr <= nm_arr))
            mn_gt = int(np.sum(mn_arr > nm_arr))
            diff = mn_arr.astype(float) - nm_arr.astype(float)
            valid = (np.std(mn_arr) > 0) and (np.std(nm_arr) > 0)
            pearson = float(np.corrcoef(mn_arr, nm_arr)[0, 1]) if valid else float("nan")

            ax.scatter(nm_arr, mn_arr, alpha=0.15, s=8, rasterized=True, edgecolors="none")
            mx = max(nm_arr.max(), mn_arr.max(), 1)
            ax.plot([0, mx], [0, mx], "k--", alpha=0.4, linewidth=1)
            ax.set_xlabel("Bhmem NM:i tag")
            ax.set_ylabel("min(4 trials)")
            ax.set_title(label, fontsize=13, fontweight="bold")
            ax.set_aspect("equal")
            ax.set_xlim(-0.5, min(mx, 50) + 0.5)
            ax.set_ylim(-0.5, min(mx, 50) + 0.5)

            stats_text = (
                f"n = {n:,}\n"
                f"Exact match: {100*exact/n:.1f}%\n"
                f"min <= NM: {100*mn_le/n:.1f}%\n"
                f"min > NM: {100*mn_gt/n:.1f}%\n"
                f"mean |diff|: {np.mean(np.abs(diff)):.2f}\n"
                f"Pearson r: {pearson:.4f}"
            )
            ax.text(0.02, -0.22, stats_text, transform=ax.transAxes,
                    fontsize=9, verticalalignment="top", fontfamily="monospace",
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.8))

        fig.suptitle("Section 1: min(converted-genome trials) vs Bhmem NM:i", fontsize=13)
        fig.subplots_adjust(bottom=0.28)
        pdf_path = f"{out_prefix}.section1_scatter.pdf"
        fig.savefig(pdf_path, dpi=150)
        plt.close(fig)
        print(f"Wrote {pdf_path}")


def section2_cross_pipeline(bhmem_results, yap_results, out_prefix):
    """Section 2: min(trials) bhmem vs min(trials) yap on shared reads."""
    shared_keys = set(bhmem_results.keys()) & set(yap_results.keys())
    if not shared_keys:
        print("Section 2: no shared reads found between bhmem and yap.")
        return

    r1_data = []
    r2_data = []
    for key in shared_keys:
        bv = bhmem_results[key]
        yv = yap_results[key]
        is_r2 = key[1]
        entry = (bv["min_trial"], yv["min_trial"], bv["mapq"], yv["mapq"],
                 bv["nm_tag"], yv["nm_tag"])
        if is_r2:
            r2_data.append(entry)
        else:
            r1_data.append(entry)

    lines = []
    lines.append("=" * 60)
    lines.append("Section 2: min(trials) bhmem vs yap (shared reads)")
    lines.append("=" * 60)
    lines.append(f"Shared reads: {len(shared_keys)}")
    lines.append("")

    for label, data in [("R1", r1_data), ("R2", r2_data), ("All", r1_data + r2_data)]:
        if not data:
            continue
        b_arr = np.array([d[0] for d in data], dtype=float)
        y_arr = np.array([d[1] for d in data], dtype=float)
        n = len(data)

        exact = int(np.sum(b_arr == y_arr))
        b_lower = int(np.sum(b_arr < y_arr))
        y_lower = int(np.sum(y_arr < b_arr))
        diff = b_arr - y_arr

        lines.append(f"--- {label} ({n} reads) ---")
        lines.append(f"  Exact match:            {exact:7d}  ({100*exact/n:.2f}%)")
        lines.append(f"  bhmem lower (better):   {b_lower:7d}  ({100*b_lower/n:.2f}%)")
        lines.append(f"  yap lower (better):     {y_lower:7d}  ({100*y_lower/n:.2f}%)")
        lines.append(f"  mean(bhmem - yap):      {np.mean(diff):+.3f}")
        lines.append(f"  median(bhmem - yap):    {np.median(diff):+.1f}")
        lines.append(f"  mean |bhmem - yap|:     {np.mean(np.abs(diff)):.3f}")
        lines.append(f"  mean bhmem min_trial:   {np.mean(b_arr):.2f}")
        lines.append(f"  mean yap min_trial:     {np.mean(y_arr):.2f}")
        if n > 1 and np.std(b_arr) > 0 and np.std(y_arr) > 0:
            pearson = float(np.corrcoef(b_arr, y_arr)[0, 1])
            lines.append(f"  Pearson r:              {pearson:.4f}")
        lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    summary_path = f"{out_prefix}.section2_summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary_text + "\n")
    print(f"Wrote {summary_path}")

    # TSV
    tsv_path = f"{out_prefix}.section2.per_read.tsv"
    with open(tsv_path, "w") as f:
        f.write("read_name\tis_r2\tmin_trial_bhmem\tmin_trial_yap\t"
                "mapq_bhmem\tmapq_yap\tnm_tag_bhmem\tnm_tag_yap\tdiff\n")
        for key in sorted(shared_keys):
            bv = bhmem_results[key]
            yv = yap_results[key]
            diff = bv["min_trial"] - yv["min_trial"]
            nm_b = bv["nm_tag"] if bv["nm_tag"] is not None else ""
            nm_y = yv["nm_tag"] if yv["nm_tag"] is not None else ""
            f.write(f"{key[0]}\t{int(key[1])}\t{bv['min_trial']}\t{yv['min_trial']}\t"
                    f"{bv['mapq']}\t{yv['mapq']}\t{nm_b}\t{nm_y}\t{diff}\n")
    print(f"Wrote {tsv_path}")

    # Scatter plots
    if HAS_MPL and (r1_data or r2_data):
        fig, axes = plt.subplots(1, 3, figsize=(18, 8))
        for ax, label, data in zip(axes, ["R1", "R2", "All"],
                                    [r1_data, r2_data, r1_data + r2_data]):
            if not data:
                ax.set_title(f"{label}: no data")
                continue
            b_arr = np.array([d[0] for d in data], dtype=float)
            y_arr = np.array([d[1] for d in data], dtype=float)
            n = len(data)
            exact = int(np.sum(b_arr == y_arr))
            b_lower = int(np.sum(b_arr < y_arr))
            y_lower = int(np.sum(y_arr < b_arr))
            diff = b_arr - y_arr
            valid = (np.std(b_arr) > 0) and (np.std(y_arr) > 0)
            pearson = float(np.corrcoef(b_arr, y_arr)[0, 1]) if valid else float("nan")

            ax.scatter(y_arr, b_arr, alpha=0.15, s=8, rasterized=True, edgecolors="none")
            mx = max(b_arr.max(), y_arr.max(), 1)
            ax.plot([0, mx], [0, mx], "k--", alpha=0.4, linewidth=1)
            ax.set_xlabel("yap min(trials)")
            ax.set_ylabel("bhmem min(trials)")
            ax.set_title(label, fontsize=13, fontweight="bold")
            ax.set_aspect("equal")
            ax.set_xlim(-0.5, min(mx, 50) + 0.5)
            ax.set_ylim(-0.5, min(mx, 50) + 0.5)

            stats_text = (
                f"n = {n:,}\n"
                f"Exact match: {100*exact/n:.1f}%\n"
                f"bhmem lower (better): {100*b_lower/n:.1f}%\n"
                f"yap lower (better): {100*y_lower/n:.1f}%\n"
                f"mean bhmem: {np.mean(b_arr):.2f}  mean yap: {np.mean(y_arr):.2f}\n"
                f"mean(bhmem-yap): {np.mean(diff):+.3f}\n"
                f"mean |diff|: {np.mean(np.abs(diff)):.2f}\n"
                f"Pearson r: {pearson:.4f}"
            )
            ax.text(0.02, -0.22, stats_text, transform=ax.transAxes,
                    fontsize=9, verticalalignment="top", fontfamily="monospace",
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="lightcyan", alpha=0.8))

        fig.suptitle("Section 2: min(trials) bhmem vs yap (shared reads)", fontsize=13)
        fig.subplots_adjust(bottom=0.32)
        pdf_path = f"{out_prefix}.section2_scatter.pdf"
        fig.savefig(pdf_path, dpi=150)
        plt.close(fig)
        print(f"Wrote {pdf_path}")


def _plot_scatter_row(axes_row, datasets, xlabel, ylabel, box_color, extra_stats_fn=None):
    """Helper to plot a row of 3 scatter panels (R1, R2, All) with stat boxes."""
    for col, (label, data) in enumerate(datasets):
        ax = axes_row[col]
        if not data:
            ax.set_title(f"{label}: no data")
            continue
        x_arr = np.array([d[0] for d in data], dtype=float)
        y_arr = np.array([d[1] for d in data], dtype=float)
        n = len(data)
        exact = int(np.sum(x_arr == y_arr))
        x_lower = int(np.sum(x_arr < y_arr))
        y_lower = int(np.sum(y_arr < x_arr))
        diff = y_arr - x_arr
        valid = n > 1 and np.std(x_arr) > 0 and np.std(y_arr) > 0
        pearson = float(np.corrcoef(x_arr, y_arr)[0, 1]) if valid else float("nan")

        ax.scatter(x_arr, y_arr, alpha=0.15, s=8, rasterized=True, edgecolors="none")
        mx = max(x_arr.max(), y_arr.max(), 1)
        ax.plot([0, mx], [0, mx], "k--", alpha=0.4, linewidth=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(label, fontsize=13, fontweight="bold")
        ax.set_aspect("equal")
        ax.set_xlim(-0.5, min(mx, 50) + 0.5)
        ax.set_ylim(-0.5, min(mx, 50) + 0.5)

        if extra_stats_fn is not None:
            stats_text = extra_stats_fn(data, n, exact, x_lower, y_lower, diff, pearson)
        else:
            stats_text = (
                f"n = {n:,}\n"
                f"Exact match: {100*exact/n:.1f}%\n"
                f"Pearson r: {pearson:.4f}"
            )
        ax.text(0.02, -0.18, stats_text, transform=ax.transAxes,
                fontsize=9, verticalalignment="top", fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor=box_color, alpha=0.8))


def combined_figure(bhmem_results, yap_results, out_prefix):
    """Combined 3x3 figure: Section 1 (row 1) + Section 2 with placement
    breakdown (row 2) + Section 3 discrepant MAPQ bhmem<30 & yap>=30 (row 3)."""
    if not HAS_MPL:
        return

    # Prepare Section 1 data
    s1_r1, s1_r2 = [], []
    for (name, is_r2), v in bhmem_results.items():
        if v["nm_tag"] is None:
            continue
        # (x=NM tag, y=min_trial)
        entry = (v["nm_tag"], v["min_trial"])
        (s1_r2 if is_r2 else s1_r1).append(entry)

    # Prepare Section 2 + 3 data
    shared_keys = set(bhmem_results.keys()) & set(yap_results.keys())
    s2_r1, s2_r2 = [], []
    s3_r1, s3_r2 = [], []  # discrepant MAPQ: bhmem<30, yap>=30
    for key in shared_keys:
        bv = bhmem_results[key]
        yv = yap_results[key]
        same_pos = (bv["chrom"] == yv["chrom"] and bv["pos"] == yv["pos"])
        same_cigar = same_pos and (bv["cigar"] == yv["cigar"])
        # (x=yap min_trial, y=bhmem min_trial, same_pos, same_cigar)
        entry = (yv["min_trial"], bv["min_trial"], same_pos, same_cigar)
        (s2_r2 if key[1] else s2_r1).append(entry)

        # Section 3: bhmem MAPQ < 30 and yap MAPQ >= 30
        if bv["mapq"] < 30 and yv["mapq"] >= 30:
            (s3_r2 if key[1] else s3_r1).append(entry)

    fig, axes = plt.subplots(3, 3, figsize=(18, 21))

    # ---- Row 1: Section 1 (min trials vs NM tag) ----
    def s1_stats(data, n, exact, x_lower, y_lower, diff, pearson):
        mn_le = int(np.sum(np.array([d[1] for d in data]) <= np.array([d[0] for d in data])))
        mn_gt = n - mn_le
        return (
            f"n = {n:,}\n"
            f"Exact match (min==NM): {100*exact/n:.1f}%\n"
            f"min <= NM: {100*mn_le/n:.1f}%\n"
            f"min > NM: {100*mn_gt/n:.1f}%\n"
            f"mean |diff|: {np.mean(np.abs(diff)):.2f}\n"
            f"Pearson r: {pearson:.4f}"
        )

    _plot_scatter_row(
        axes[0],
        [("R1", s1_r1), ("R2", s1_r2), ("All", s1_r1 + s1_r2)],
        xlabel="Bhmem NM:i tag", ylabel="min(4 trials)",
        box_color="wheat", extra_stats_fn=s1_stats,
    )

    # ---- Row 2: Section 2 (bhmem vs yap, all shared reads) ----
    def s2_stats(data, n, exact, x_lower, y_lower, diff, pearson):
        b_arr = np.array([d[1] for d in data], dtype=float)
        y_arr = np.array([d[0] for d in data], dtype=float)
        b_lower = int(np.sum(b_arr < y_arr))
        y_lower_ct = int(np.sum(y_arr < b_arr))
        same_pos_ct = sum(1 for d in data if d[2])
        same_cigar_ct = sum(1 for d in data if d[3])
        return (
            f"n = {n:,}\n"
            f"Same min(trials): {100*exact/n:.1f}%\n"
            f"  same pos+cigar: {100*same_cigar_ct/n:.1f}%\n"
            f"  same pos only:  {100*same_pos_ct/n:.1f}%\n"
            f"bhmem lower: {100*b_lower/n:.1f}%\n"
            f"yap lower: {100*y_lower_ct/n:.1f}%\n"
            f"mean bhmem: {np.mean(b_arr):.2f}  yap: {np.mean(y_arr):.2f}\n"
            f"Pearson r: {pearson:.4f}"
        )

    _plot_scatter_row(
        axes[1],
        [("R1", s2_r1), ("R2", s2_r2), ("All", s2_r1 + s2_r2)],
        xlabel="yap min(trials)", ylabel="bhmem min(trials)",
        box_color="lightcyan", extra_stats_fn=s2_stats,
    )

    # ---- Row 3: Section 3 (discrepant MAPQ: bhmem<30, yap>=30) ----
    def s3_stats(data, n, exact, x_lower, y_lower, diff, pearson):
        b_arr = np.array([d[1] for d in data], dtype=float)
        y_arr = np.array([d[0] for d in data], dtype=float)
        b_lower = int(np.sum(b_arr < y_arr))
        y_lower_ct = int(np.sum(y_arr < b_arr))
        same_pos_ct = sum(1 for d in data if d[2])
        same_cigar_ct = sum(1 for d in data if d[3])
        return (
            f"n = {n:,}\n"
            f"Filter: bhmem MQ<30, yap MQ>=30\n"
            f"Same min(trials): {100*exact/n:.1f}%\n"
            f"  same pos+cigar: {100*same_cigar_ct/n:.1f}%\n"
            f"  same pos only:  {100*same_pos_ct/n:.1f}%\n"
            f"bhmem lower: {100*b_lower/n:.1f}%\n"
            f"yap lower: {100*y_lower_ct/n:.1f}%\n"
            f"mean bhmem: {np.mean(b_arr):.2f}  yap: {np.mean(y_arr):.2f}\n"
            f"Pearson r: {pearson:.4f}"
        )

    _plot_scatter_row(
        axes[2],
        [("R1", s3_r1), ("R2", s3_r2), ("All", s3_r1 + s3_r2)],
        xlabel="yap min(trials)", ylabel="bhmem min(trials)",
        box_color="mistyrose", extra_stats_fn=s3_stats,
    )

    # Row labels
    axes[0, 0].annotate("Section 1\nmin(trials) vs\nBhmem NM:i",
                        xy=(-0.35, 0.5), xycoords="axes fraction",
                        fontsize=11, fontweight="bold", ha="center", va="center",
                        rotation=90)
    axes[1, 0].annotate("Section 2\nbhmem vs yap\n(all shared)",
                        xy=(-0.35, 0.5), xycoords="axes fraction",
                        fontsize=11, fontweight="bold", ha="center", va="center",
                        rotation=90)
    axes[2, 0].annotate("Section 3\nbhmem MQ<30\nyap MQ>=30",
                        xy=(-0.35, 0.5), xycoords="axes fraction",
                        fontsize=11, fontweight="bold", ha="center", va="center",
                        rotation=90)

    fig.subplots_adjust(left=0.12, bottom=0.08, hspace=0.55, wspace=0.3)
    pdf_path = f"{out_prefix}.combined.pdf"
    fig.savefig(pdf_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {pdf_path}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("bhmem_bam", help="Bhmem BAM file")
    ap.add_argument("yap_bam", help="Yap 3C sorted BAM file")
    ap.add_argument(
        "bisulfite_genome",
        nargs="?",
        default=DEFAULT_BISULFITE,
        help=f"Bismark Bisulfite_Genome dir (default: {DEFAULT_BISULFITE})",
    )
    ap.add_argument(
        "--joined-tsv",
        help="mapq_comparison joined TSV to restrict to shared reads",
    )
    ap.add_argument("-o", "--output-prefix", default="min_trial_nm_report",
                    help="Output prefix for TSV/PDF/summary")
    ap.add_argument("--max-reads", type=int, default=20000,
                    help="Max primary reads per BAM (0 = all)")
    ap.add_argument("--max-joined-rows", type=int, default=0,
                    help="Max rows from joined TSV (0 = all)")
    args = ap.parse_args()

    ct_path = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_path = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for p in (ct_path, ga_path, args.bhmem_bam, args.yap_bam):
        if not os.path.isfile(p):
            print(f"ERROR: missing file: {p}", file=sys.stderr)
            sys.exit(1)

    ct_fa = pysam.FastaFile(ct_path)
    ga_fa = pysam.FastaFile(ga_path)

    # Step 1: Scan yap BAM to get all yap read keys (fast — no trial computation)
    print(f"\n--- Scanning yap BAM for read keys ---", flush=True)
    yap_all_keys = set()
    yap_bam_scan = pysam.AlignmentFile(args.yap_bam, "rb")
    for read in yap_bam_scan:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        parsed = _normalize_yap_key(read.query_name)
        if parsed is not None:
            yap_all_keys.add(parsed)
    yap_bam_scan.close()
    print(f"  {len(yap_all_keys)} unique yap read keys", flush=True)

    # Step 2: Compute min(trials) on bhmem BAM — only shared reads
    print(f"\n--- Computing min(trials) on bhmem BAM (shared reads only) ---", flush=True)
    bhmem_results = compute_min_trials(
        args.bhmem_bam, ct_fa, ga_fa,
        read_set=yap_all_keys, max_reads=args.max_reads,
    )
    print(f"  {len(bhmem_results)} reads processed", flush=True)

    # Section 1: bhmem min(trials) vs NM tag
    print(f"\n", flush=True)
    section1_bhmem_agreement(bhmem_results, args.output_prefix)

    # Step 3: Compute min(trials) on yap BAM — restrict to bhmem result keys
    bhmem_keys = set(bhmem_results.keys())
    print(f"\n--- Computing min(trials) on yap BAM (restricting to {len(bhmem_keys)} bhmem keys) ---", flush=True)
    yap_results = compute_min_trials(
        args.yap_bam, ct_fa, ga_fa,
        read_set=bhmem_keys, max_reads=0,  # no limit — find all matches
        name_mode="yap",
    )
    print(f"  {len(yap_results)} reads processed", flush=True)

    # Section 2: cross-pipeline
    print(f"\n", flush=True)
    section2_cross_pipeline(bhmem_results, yap_results, args.output_prefix)

    # Combined figure
    print(f"\n", flush=True)
    combined_figure(bhmem_results, yap_results, args.output_prefix)

    ct_fa.close()
    ga_fa.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
