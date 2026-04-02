#!/usr/bin/env python3
"""
Cross-pipeline bisulfite-aware NM comparison: bhmem NM:i tag vs yap recomputed NM.

- **Bhmem**: use ``NM:i`` tag directly (already bisulfite-aware from BWA converted index).
- **Yap (Bismark)**: use ``XR:Z`` (read conversion) and ``XG:Z`` (genome conversion) tags
  to select the correct converted reference + query conversion, walk CIGAR to compute NM.

Both NM values are in the same space: edit distance of converted query vs converted reference.

Outputs: per-read TSV, summary stats, scatter PDF (R1/R2/All).

Example::

  python report_bisulfite_nm_cross_pipeline.py \\
    /path/to/bhmem.bam \\
    /path/to/yap.3C.sorted.bam \\
    /path/to/Bisulfite_Genome \\
    -o output_prefix --max-reads 20000
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import (
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
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


def _normalize_yap_key(qname):
    """Parse yap QNAME -> (base_id, is_r2) or None."""
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


def _apply_conversion(seq, conversion):
    """Apply bisulfite conversion to a sequence.

    conversion: 'CT' -> C->T, 'GA' -> G->A
    """
    s = seq.upper()
    if conversion == "CT":
        return s.replace("C", "T")
    elif conversion == "GA":
        return s.replace("G", "A")
    return s


def compute_yap_bisulfite_nm(read, ct_fa, ga_fa):
    """Compute bisulfite-aware NM for a yap (Bismark) read using XR/XG tags.

    XR = read conversion applied by Bismark (CT or GA)
    XG = genome strand (CT or GA genome)

    The read was aligned as: convert(query) vs converted_genome.
    We reconstruct that and walk the CIGAR.

    Returns (nm, method) or (-1, reason).
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0, "unmapped"

    if not read.has_tag("XR") or not read.has_tag("XG"):
        return -1, "no_XR_XG"

    xr = str(read.get_tag("XR")).upper()  # read conversion: CT or GA
    xg = str(read.get_tag("XG")).upper()  # genome: CT or GA

    if xr not in ("CT", "GA") or xg not in ("CT", "GA"):
        return -1, "bad_XR_XG"

    # Select converted reference
    fa = ct_fa if xg == "CT" else ga_fa
    ref_contig = bisulfite_converted_contig_name(fa, read.reference_name, xg)
    if ref_contig is None:
        return -1, "no_contig"

    # Build converted query
    seq = read.query_sequence.upper()
    query_conv = _apply_conversion(seq, xr)

    # Walk CIGAR
    nm = count_nm_style_edit_distance_converted_explicit(
        read, fa, ref_contig=ref_contig, query_converted=query_conv,
    )
    if nm < 0:
        # Try RC variant (Bismark may store SEQ in different orientation)
        query_conv_v2 = _reverse_complement_dna(_apply_conversion(_reverse_complement_dna(seq), xr))
        nm2 = count_nm_style_edit_distance_converted_explicit(
            read, fa, ref_contig=ref_contig, query_converted=query_conv_v2,
        )
        if nm2 >= 0:
            return nm2, "xr_xg_v2"
        return -1, "cigar_fail"

    return nm, "xr_xg"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("bhmem_bam", help="Bhmem BAM file")
    ap.add_argument("yap_bam", help="Yap 3C sorted BAM file")
    ap.add_argument(
        "bisulfite_genome", nargs="?", default=DEFAULT_BISULFITE,
        help=f"Bismark Bisulfite_Genome dir (default: {DEFAULT_BISULFITE})",
    )
    ap.add_argument("-o", "--output-prefix", default="bisulfite_nm_cross")
    ap.add_argument("--max-reads", type=int, default=20000,
                    help="Max shared reads to process (0 = all)")
    args = ap.parse_args()

    ct_path = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_path = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for p in (ct_path, ga_path, args.bhmem_bam, args.yap_bam):
        if not os.path.isfile(p):
            print(f"ERROR: missing: {p}", file=sys.stderr)
            sys.exit(1)

    ct_fa = pysam.FastaFile(ct_path)
    ga_fa = pysam.FastaFile(ga_path)

    # Step 1: Scan bhmem BAM — collect NM:i for all reads
    print("Step 1: Reading bhmem NM:i tags...", flush=True)
    bhmem_data = {}  # (base_id, is_r2) -> {nm, mapq, chrom, pos}
    bam_b = pysam.AlignmentFile(args.bhmem_bam, "rb")
    for read in bam_b:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("NM"):
            continue
        is_r2 = bool(read.is_paired and read.is_read2)
        key = (read.query_name, is_r2)
        bhmem_data[key] = {
            "nm": int(read.get_tag("NM")),
            "mapq": int(read.mapping_quality),
            "chrom": read.reference_name,
            "pos": read.reference_start,
        }
    bam_b.close()
    print(f"  {len(bhmem_data)} bhmem reads", flush=True)

    # Step 2: Scan yap BAM — compute bisulfite NM for shared reads
    print("Step 2: Computing yap bisulfite NM (shared reads)...", flush=True)
    yap_data = {}  # (base_id, is_r2) -> {nm, mapq, chrom, pos, method}
    bam_y = pysam.AlignmentFile(args.yap_bam, "rb")
    n_processed = 0
    n_shared = 0
    method_counts = {}

    for read in bam_y:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None:
            continue

        parsed = _normalize_yap_key(read.query_name)
        if parsed is None:
            continue
        key = parsed

        if key not in bhmem_data:
            continue

        # For split reads, keep the one with lowest NM
        nm, method = compute_yap_bisulfite_nm(read, ct_fa, ga_fa)
        mapq = int(read.mapping_quality)

        n_processed += 1
        method_counts[method] = method_counts.get(method, 0) + 1

        if nm < 0:
            continue

        entry = {
            "nm": nm,
            "mapq": mapq,
            "chrom": read.reference_name,
            "pos": read.reference_start,
            "method": method,
        }

        if key in yap_data:
            if nm < yap_data[key]["nm"]:
                yap_data[key] = entry
        else:
            yap_data[key] = entry
            n_shared += 1

        if args.max_reads and n_shared >= args.max_reads:
            break

        if n_shared % 5000 == 0 and n_shared > 0:
            print(f"  {n_shared} shared reads...", flush=True)

    bam_y.close()
    print(f"  {n_processed} yap reads scanned, {n_shared} shared with valid NM", flush=True)
    print(f"  Methods: {method_counts}", flush=True)

    ct_fa.close()
    ga_fa.close()

    # Step 3: Join and compare
    shared_keys = set(bhmem_data.keys()) & set(yap_data.keys())
    print(f"\nJoined: {len(shared_keys)} shared reads\n", flush=True)

    if not shared_keys:
        print("No shared reads found.")
        return

    r1_data = []  # (bhmem_nm, yap_nm, bhmem_mq, yap_mq, same_pos)
    r2_data = []

    for key in shared_keys:
        bv = bhmem_data[key]
        yv = yap_data[key]
        same_pos = (bv["chrom"] == yv["chrom"] and abs(bv["pos"] - yv["pos"]) <= 5)
        entry = (bv["nm"], yv["nm"], bv["mapq"], yv["mapq"], same_pos)
        if key[1]:
            r2_data.append(entry)
        else:
            r1_data.append(entry)

    # Print and save summary
    lines = []
    lines.append("=" * 60)
    lines.append("Bisulfite-aware NM: bhmem tag vs yap XR/XG recompute")
    lines.append("=" * 60)
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
        same_pos_ct = sum(1 for d in data if d[4])
        valid = n > 1 and np.std(b_arr) > 0 and np.std(y_arr) > 0
        pearson = float(np.corrcoef(b_arr, y_arr)[0, 1]) if valid else float("nan")

        # Among same-position reads
        sp_data = [d for d in data if d[4]]
        sp_n = len(sp_data)
        if sp_n > 0:
            sp_b = np.array([d[0] for d in sp_data], dtype=float)
            sp_y = np.array([d[1] for d in sp_data], dtype=float)
            sp_exact = int(np.sum(sp_b == sp_y))
        else:
            sp_exact = 0

        # Among reads with bhmem MQ<30, yap MQ>=30
        disc_data = [d for d in data if d[2] < 30 and d[3] >= 30]
        disc_n = len(disc_data)
        if disc_n > 0:
            disc_b = np.array([d[0] for d in disc_data], dtype=float)
            disc_y = np.array([d[1] for d in disc_data], dtype=float)
            disc_b_lower = int(np.sum(disc_b < disc_y))
            disc_y_lower = int(np.sum(disc_y < disc_b))
        else:
            disc_b_lower = disc_y_lower = 0

        lines.append(f"--- {label} ({n} reads) ---")
        lines.append(f"  Same NM:                {exact:7d}  ({100*exact/n:.2f}%)")
        lines.append(f"  bhmem NM lower (better):{b_lower:7d}  ({100*b_lower/n:.2f}%)")
        lines.append(f"  yap NM lower (better):  {y_lower:7d}  ({100*y_lower/n:.2f}%)")
        lines.append(f"  mean bhmem NM:          {np.mean(b_arr):.2f}")
        lines.append(f"  mean yap NM:            {np.mean(y_arr):.2f}")
        lines.append(f"  mean(bhmem - yap):      {np.mean(diff):+.3f}")
        lines.append(f"  mean |diff|:            {np.mean(np.abs(diff)):.3f}")
        lines.append(f"  Pearson r:              {pearson:.4f}")
        lines.append(f"  Same position (±5bp):   {same_pos_ct:7d}  ({100*same_pos_ct/n:.2f}%)")
        if sp_n > 0:
            lines.append(f"    same NM at same pos:  {sp_exact:7d}  ({100*sp_exact/sp_n:.2f}%)")
        lines.append(f"  Discrepant MQ (b<30,y>=30): {disc_n}")
        if disc_n > 0:
            lines.append(f"    bhmem lower: {disc_b_lower}  yap lower: {disc_y_lower}")
        lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    summary_path = f"{args.output_prefix}.summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary_text + "\n")
    print(f"Wrote {summary_path}")

    # TSV
    tsv_path = f"{args.output_prefix}.per_read.tsv"
    with open(tsv_path, "w") as f:
        f.write("read_name\tis_r2\tbhmem_nm\tyap_nm\tbhmem_mapq\tyap_mapq\t"
                "same_pos\tbhmem_chrom\tbhmem_pos\tyap_chrom\tyap_pos\tdiff\n")
        for key in sorted(shared_keys):
            bv = bhmem_data[key]
            yv = yap_data[key]
            same_pos = (bv["chrom"] == yv["chrom"] and abs(bv["pos"] - yv["pos"]) <= 5)
            diff = bv["nm"] - yv["nm"]
            f.write(f"{key[0]}\t{int(key[1])}\t{bv['nm']}\t{yv['nm']}\t"
                    f"{bv['mapq']}\t{yv['mapq']}\t{int(same_pos)}\t"
                    f"{bv['chrom']}\t{bv['pos']}\t{yv['chrom']}\t{yv['pos']}\t{diff}\n")
    print(f"Wrote {tsv_path}")

    # Combined 3x3 scatter PDF
    if HAS_MPL:
        fig, axes = plt.subplots(3, 3, figsize=(18, 21))

        def _plot_row(row_axes, datasets, row_title, box_color, filter_label=None):
            for col, (label, data) in enumerate(datasets):
                ax = row_axes[col]
                if not data:
                    ax.set_title(f"{label}: no data")
                    continue
                b_arr = np.array([d[0] for d in data], dtype=float)
                y_arr = np.array([d[1] for d in data], dtype=float)
                n = len(data)
                exact = int(np.sum(b_arr == y_arr))
                b_lower = int(np.sum(b_arr < y_arr))
                y_lower = int(np.sum(y_arr < b_arr))
                same_pos_ct = sum(1 for d in data if d[4])
                valid = n > 1 and np.std(b_arr) > 0 and np.std(y_arr) > 0
                pearson = float(np.corrcoef(b_arr, y_arr)[0, 1]) if valid else float("nan")

                ax.scatter(y_arr, b_arr, alpha=0.15, s=8, rasterized=True, edgecolors="none")
                mx = max(b_arr.max(), y_arr.max(), 1)
                ax.plot([0, mx], [0, mx], "k--", alpha=0.4, linewidth=1)
                ax.set_xlabel("yap NM (XR/XG recompute)")
                ax.set_ylabel("bhmem NM:i tag")
                ax.set_title(label, fontsize=13, fontweight="bold")
                ax.set_aspect("equal")
                ax.set_xlim(-0.5, min(mx, 50) + 0.5)
                ax.set_ylim(-0.5, min(mx, 50) + 0.5)

                stats = (
                    f"n = {n:,}\n"
                    f"Same NM: {100*exact/n:.1f}%\n"
                    f"bhmem lower: {100*b_lower/n:.1f}%\n"
                    f"yap lower: {100*y_lower/n:.1f}%\n"
                    f"mean bhmem: {np.mean(b_arr):.2f}  yap: {np.mean(y_arr):.2f}\n"
                    f"same pos: {100*same_pos_ct/n:.1f}%\n"
                    f"Pearson r: {pearson:.4f}"
                )
                if filter_label:
                    stats = f"{filter_label}\n{stats}"
                ax.text(0.02, -0.18, stats, transform=ax.transAxes,
                        fontsize=9, verticalalignment="top", fontfamily="monospace",
                        bbox=dict(boxstyle="round,pad=0.4", facecolor=box_color, alpha=0.8))

        # Row 1: All shared reads
        _plot_row(axes[0],
                  [("R1", r1_data), ("R2", r2_data), ("All", r1_data + r2_data)],
                  "All shared reads", "lightcyan")
        axes[0, 0].annotate("All shared\nreads",
                            xy=(-0.35, 0.5), xycoords="axes fraction",
                            fontsize=11, fontweight="bold", ha="center", va="center",
                            rotation=90)

        # Row 2: Same position only
        sp_r1 = [d for d in r1_data if d[4]]
        sp_r2 = [d for d in r2_data if d[4]]
        _plot_row(axes[1],
                  [("R1", sp_r1), ("R2", sp_r2), ("All", sp_r1 + sp_r2)],
                  "Same position", "wheat")
        axes[1, 0].annotate("Same position\n(±5bp)",
                            xy=(-0.35, 0.5), xycoords="axes fraction",
                            fontsize=11, fontweight="bold", ha="center", va="center",
                            rotation=90)

        # Row 3: Discrepant MAPQ (bhmem<30, yap>=30)
        disc_r1 = [d for d in r1_data if d[2] < 30 and d[3] >= 30]
        disc_r2 = [d for d in r2_data if d[2] < 30 and d[3] >= 30]
        _plot_row(axes[2],
                  [("R1", disc_r1), ("R2", disc_r2), ("All", disc_r1 + disc_r2)],
                  "Discrepant MAPQ", "mistyrose",
                  filter_label="bhmem MQ<30, yap MQ>=30")
        axes[2, 0].annotate("bhmem MQ<30\nyap MQ>=30",
                            xy=(-0.35, 0.5), xycoords="axes fraction",
                            fontsize=11, fontweight="bold", ha="center", va="center",
                            rotation=90)

        fig.suptitle("Bisulfite-aware NM: bhmem NM:i tag vs yap XR/XG recompute",
                     fontsize=14)
        fig.subplots_adjust(left=0.12, bottom=0.06, hspace=0.55, wspace=0.3)
        pdf_path = f"{args.output_prefix}.pdf"
        fig.savefig(pdf_path, dpi=150)
        plt.close(fig)
        print(f"Wrote {pdf_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
