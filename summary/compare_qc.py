#!/usr/bin/env python
"""
Compare QC metrics between two scNOME/m3C datasets.

Reads MappingSummary.csv.gz from two directories and produces:
  1. A combined comparison CSV
  2. Side-by-side violin/box plots for key metrics
  3. A summary statistics table

Usage:
    python compare_qc.py \
        --dataset1 /path/to/methylhic/alignment/stats/MappingSummary.csv.gz \
        --dataset2 /path/to/snmCseq3/alignment/stats/MappingSummary.csv.gz \
        --label1 "MethylHiC" \
        --label2 "snmC-seq3" \
        --output compare_qc

Outputs:
    compare_qc.combined.csv
    compare_qc.summary.csv
    compare_qc.plots.pdf
"""

import argparse
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════════════
# Define which metrics to compare and how to group them
# ═══════════════════════════════════════════════════════════════════

METRIC_GROUPS = {
    "Read Counts": [
        ("CellInputReadPairs", "Input Read Pairs"),
        ("R1TrimmedReads", "R1 Trimmed Reads"),
        ("R2TrimmedReads", "R2 Trimmed Reads"),
        ("FinalmCReads", "Final mC Reads"),
    ],
    "Mapping Rate (%)": [
        ("R1MappingRate", "R1 Mapping Rate"),
        ("R2MappingRate", "R2 Mapping Rate"),
        ("R1MappingRateMapQ30", "R1 MapQ30 Rate"),
        ("R2MappingRateMapQ30", "R2 MapQ30 Rate"),
    ],
    "Duplication Rate (%)": [
        ("R1DuplicationRate", "R1 Dup Rate"),
        ("R2DuplicationRate", "R2 Dup Rate"),
    ],
    "Contacts": [
        ("TotalContacts", "Total Contacts"),
        ("CisLongContact", "Cis Long Contacts"),
        ("TransContact", "Trans Contacts"),
    ],
    "Contact Ratios": [
        ("CisShortRatio", "Cis Short Ratio"),
        ("CisLongRatio", "Cis Long Ratio"),
        ("TransRatio", "Trans Ratio"),
    ],
    "Methylation": [
        ("mCGFrac", "mCG Fraction"),
        ("mCHFrac", "mCH Fraction"),
        ("mCCCFrac", "mCCC Fraction"),
        ("LambdaCYFrac", "Lambda CY Fraction"),
    ],
    "Coverage": [
        ("GenomeCov", "Genome Coverage"),
    ],
}


# ═══════════════════════════════════════════════════════════════════
# Load and label datasets
# ═══════════════════════════════════════════════════════════════════

def load_dataset(path, label):
    df = pd.read_csv(path, index_col=0)
    df.insert(0, "Dataset", label)
    return df


# ═══════════════════════════════════════════════════════════════════
# Summary statistics
# ═══════════════════════════════════════════════════════════════════

def compute_summary(dfs, labels):
    """Compute summary stats for metrics shared across all datasets."""
    all_metrics = set()
    for group_metrics in METRIC_GROUPS.values():
        for col, _ in group_metrics:
            all_metrics.add(col)

    shared_cols = sorted(set.intersection(*(set(df.columns) for df in dfs)) & all_metrics)

    rows = []
    for col in shared_cols:
        row = {"Metric": col}
        for df, label in zip(dfs, labels):
            s = pd.to_numeric(df[col], errors="coerce").dropna()
            if len(s) == 0:
                row[f"{label}_N"] = 0
                row[f"{label}_Median"] = None
                row[f"{label}_Mean"] = None
                row[f"{label}_Std"] = None
                row[f"{label}_Min"] = None
                row[f"{label}_Max"] = None
            else:
                row[f"{label}_N"] = len(s)
                row[f"{label}_Median"] = s.median()
                row[f"{label}_Mean"] = s.mean()
                row[f"{label}_Std"] = s.std()
                row[f"{label}_Min"] = s.min()
                row[f"{label}_Max"] = s.max()
        rows.append(row)

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════

def plot_comparison(dfs, labels, pdf_path):
    """Generate grouped violin/box plots for N datasets."""
    color_palette = [
        "#4C78A8",  # blue
        "#E45756",  # red
        "#72B7B2",  # teal
        "#59A14F",  # green
        "#F28E2B",  # orange
        "#B279A2",  # purple
    ]
    colors = {label: color_palette[i % len(color_palette)] for i, label in enumerate(labels)}

    with PdfPages(pdf_path) as pdf:
        for group_name, metrics in METRIC_GROUPS.items():
            # Keep metrics shared across all datasets.
            valid_metrics = []
            for col, display in metrics:
                if all(col in df.columns for df in dfs):
                    # Only plot the metric if at least one dataset has values
                    any_nonempty = False
                    for df in dfs:
                        s = pd.to_numeric(df[col], errors="coerce").dropna()
                        if len(s) > 0:
                            any_nonempty = True
                            break
                    if any_nonempty:
                        valid_metrics.append((col, display))

            if not valid_metrics:
                continue

            n_metrics = len(valid_metrics)
            fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5))
            if n_metrics == 1:
                axes = [axes]

            fig.suptitle(group_name, fontsize=16, fontweight="bold", y=1.02)

            for ax, (col, display) in zip(axes, valid_metrics):
                # Build per-dataset distributions (skip empty ones).
                data_to_plot = []
                positions = []
                clrs = []
                labels_plot = []

                for i, (df, label) in enumerate(zip(dfs, labels), start=1):
                    s = pd.to_numeric(df[col], errors="coerce").dropna().values
                    if len(s) == 0:
                        continue
                    data_to_plot.append(s)
                    positions.append(len(positions) + 1)  # compact positions for non-empty datasets
                    clrs.append(colors[label])
                    labels_plot.append(label)

                # Violin plot + overlays
                if len(data_to_plot) > 0:
                    parts = ax.violinplot(
                        data_to_plot,
                        positions=positions,
                        showmeans=False,
                        showmedians=False,
                        showextrema=False,
                    )
                    for pc, color in zip(parts['bodies'], clrs):
                        pc.set_facecolor(color)
                        pc.set_alpha(0.3)

                    bp = ax.boxplot(
                        data_to_plot,
                        positions=positions,
                        widths=0.15,
                        patch_artist=True,
                        showfliers=False,
                        zorder=3,
                    )
                    for patch, color in zip(bp['boxes'], clrs):
                        patch.set_facecolor(color)
                        patch.set_alpha(0.8)
                    for element in ['whiskers', 'caps', 'medians']:
                        for line in bp[element]:
                            line.set_color('#333333')
                            line.set_linewidth(1.2)

                    for pos, vals, color in zip(positions, data_to_plot, clrs):
                        jitter = np.random.normal(0, 0.03, size=len(vals))
                        ax.scatter(
                            pos + jitter,
                            vals,
                            c=color,
                            alpha=0.4,
                            s=8,
                            zorder=2,
                            edgecolors="none",
                        )

                ax.set_title(display, fontsize=11)
                ax.set_xticks(positions)
                ax.set_xticklabels(labels_plot, fontsize=9, rotation=15)
                ax.tick_params(axis='y', labelsize=9)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                # Add median annotation
                for pos, vals in zip(positions, data_to_plot):
                    med = np.median(vals)
                    ax.annotate(
                        f"{med:.4g}",
                        xy=(pos, med),
                        xytext=(pos + 0.3, med),
                        fontsize=7,
                        color="#555555",
                        arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.5),
                    )

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ── Summary table page ──
        fig, ax = plt.subplots(figsize=(14, 8))
        ax.axis('off')
        title = "Summary: " + " vs ".join(labels)
        ax.set_title(title, fontsize=14, fontweight="bold", pad=20)

        summary_rows = []
        for group_name, metrics in METRIC_GROUPS.items():
            for col, display in metrics:
                if all(col in df.columns for df in dfs):
                    med_cells = []
                    any_nonempty = False
                    for df in dfs:
                        s = pd.to_numeric(df[col], errors="coerce").dropna()
                        if len(s) > 0:
                            any_nonempty = True
                            med_cells.append((f"{s.median():.4g}", f"{s.mean():.4g}"))
                        else:
                            med_cells.append(("N/A", "N/A"))
                    if any_nonempty:
                        row = [display]
                        for med, mean in med_cells:
                            row.extend([med, mean])
                        summary_rows.append(row)

        if summary_rows:
            colLabels = ["Metric"]
            for label in labels:
                colLabels.extend([f"{label}\nMedian", f"{label}\nMean"])

            table = ax.table(
                cellText=summary_rows,
                colLabels=colLabels,
                loc="center",
                cellLoc="center",
            )
            table.auto_set_font_size(False)
            table.set_fontsize(7)
            table.scale(1.0, 1.3)

            # Style header
            for (row, col), cell in table.get_celld().items():
                if row == 0:
                    cell.set_facecolor("#4C78A8")
                    cell.set_text_props(color="white", fontweight="bold")
                elif row % 2 == 0:
                    cell.set_facecolor("#f0f4f8")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"Saved plots: {pdf_path}")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Compare QC metrics between two scNOME/m3C datasets"
    )
    parser.add_argument("--dataset1", required=True,
                        help="Path to first MappingSummary.csv.gz")
    parser.add_argument("--dataset2", required=True,
                        help="Path to second MappingSummary.csv.gz")
    parser.add_argument("--label1", default="Dataset1",
                        help="Label for first dataset")
    parser.add_argument("--label2", default="Dataset2",
                        help="Label for second dataset")
    parser.add_argument("--dataset3", required=False, default=None,
                        help="Optional Path to third MappingSummary.csv.gz")
    parser.add_argument("--dataset4", required=False, default=None,
                        help="Optional Path to fourth MappingSummary.csv.gz")
    parser.add_argument("--label3", default="Dataset3",
                        help="Label for third dataset (used if --dataset3 is set)")
    parser.add_argument("--label4", default="Dataset4",
                        help="Label for fourth dataset (used if --dataset4 is set)")
    parser.add_argument("--output", default="compare_qc",
                        help="Output prefix")
    args = parser.parse_args()

    datasets = [(args.dataset1, args.label1), (args.dataset2, args.label2)]
    if args.dataset3 is not None:
        datasets.append((args.dataset3, args.label3))
    if args.dataset4 is not None:
        datasets.append((args.dataset4, args.label4))

    dfs = []
    labels = []
    for path, label in datasets:
        print(f"Dataset: {path} ({label})")
        df = load_dataset(path, label)
        dfs.append(df)
        labels.append(label)

    for df, label in zip(dfs, labels):
        print(f"  {label}: {len(df)} cells")

    # Combined CSV
    combined = pd.concat(dfs)
    combined_path = f"{args.output}.combined.csv"
    combined.to_csv(combined_path)
    print(f"\nSaved combined: {combined_path}")

    # Summary stats
    summary = compute_summary(dfs, labels)
    summary_path = f"{args.output}.summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"Saved summary: {summary_path}")

    # Print a compact summary line (median only)
    header = ["Metric".ljust(35)] + [f"{lbl} (med)".rjust(14) for lbl in labels]
    print("\n" + "  ".join(header))
    print("=" * (20 * len(labels) + 35))
    for _, row in summary.iterrows():
        parts = [str(row["Metric"]).ljust(35)]
        for lbl in labels:
            m = row.get(f"{lbl}_Median")
            parts.append(f"{m:.4g}" if pd.notna(m) else "N/A".rjust(14))
        print("  ".join(parts))

    # Plots (all datasets together)
    pdf_path = f"{args.output}.plots.pdf"
    plot_comparison(dfs, labels, pdf_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
