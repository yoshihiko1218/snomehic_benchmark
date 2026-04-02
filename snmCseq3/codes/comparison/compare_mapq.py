#!/usr/bin/env python3
"""
Compare MAPQ from bhmem (BwaMem) and yap (Bowtie2) for the same reads.

Input: two TSV files from extract_mapq.py (base_id, is_r1, mapq)
Output: joined table, summary stats, correlation, and optional plots.
"""

import argparse
import sys
from collections import defaultdict

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def load_tsv(path: str):
    """Return dict: (base_id, is_r1) -> mapq"""
    d = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            base_id, is_r1, mapq = parts[0], int(parts[1]), int(parts[2])
            d[(base_id, is_r1)] = mapq
    return d


def main():
    ap = argparse.ArgumentParser(description="Compare MAPQ between bhmem and yap")
    ap.add_argument("bhmem_tsv", help="TSV from extract_mapq.py for bhmem BAM")
    ap.add_argument("yap_tsv", help="TSV from extract_mapq.py for yap 3C BAM")
    ap.add_argument("-o", "--output-prefix", default="mapq_comparison", help="Output prefix for TSV and plots")
    ap.add_argument("--plot", action="store_true", help="Generate scatter and distribution plots")
    args = ap.parse_args()

    bhmem = load_tsv(args.bhmem_tsv)
    yap = load_tsv(args.yap_tsv)

    # Join on (base_id, is_r1)
    keys = sorted(set(bhmem) & set(yap))
    if not keys:
        print("ERROR: No overlapping reads between bhmem and yap TSVs.", file=sys.stderr)
        sys.exit(1)

    mapq_b = np.array([bhmem[k] for k in keys])
    mapq_y = np.array([yap[k] for k in keys])

    # Write joined table
    out_tsv = f"{args.output_prefix}.joined.tsv"
    with open(out_tsv, "w") as f:
        f.write("base_id\tis_r1\tmapq_bhmem\tmapq_yap\tdiff\n")
        for k, mb, my in zip(keys, mapq_b, mapq_y):
            f.write(f"{k[0]}\t{k[1]}\t{mb}\t{my}\t{mb - my}\n")
    print(f"Wrote {len(keys)} joined rows to {out_tsv}")

    # Summary stats
    diff = mapq_b - mapq_y
    stats = [
        ("N_matched", len(keys)),
        ("bhmem_mean", float(np.mean(mapq_b))),
        ("bhmem_median", float(np.median(mapq_b))),
        ("bhmem_std", float(np.std(mapq_b))),
        ("yap_mean", float(np.mean(mapq_y))),
        ("yap_median", float(np.median(mapq_y))),
        ("yap_std", float(np.std(mapq_y))),
        ("diff_mean", float(np.mean(diff))),
        ("diff_median", float(np.median(diff))),
        ("diff_std", float(np.std(diff))),
        ("correlation", float(np.corrcoef(mapq_b, mapq_y)[0, 1])),
        ("n_bhmem_higher", int(np.sum(diff > 0))),
        ("n_yap_higher", int(np.sum(diff < 0))),
        ("n_equal", int(np.sum(diff == 0))),
    ]
    stats_path = f"{args.output_prefix}.stats.txt"
    with open(stats_path, "w") as f:
        for name, val in stats:
            f.write(f"{name}\t{val}\n")
    print(f"Wrote stats to {stats_path}")
    for name, val in stats:
        print(f"  {name}: {val}")

    if args.plot and HAS_MATPLOTLIB:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))

        # Scatter
        ax = axes[0]
        ax.scatter(mapq_y, mapq_b, alpha=0.3, s=2)
        mx = max(mapq_b.max(), mapq_y.max())
        ax.plot([0, mx], [0, mx], "k--", alpha=0.5, label="y=x")
        ax.set_xlabel("yap (Bowtie2) MAPQ")
        ax.set_ylabel("bhmem (BwaMem) MAPQ")
        ax.set_title("MAPQ comparison (same reads)")
        ax.legend()
        ax.set_aspect("equal")

        # Diff histogram
        ax = axes[1]
        ax.hist(diff, bins=min(51, len(np.unique(diff)) + 1), edgecolor="black", alpha=0.7)
        ax.axvline(0, color="red", linestyle="--")
        ax.set_xlabel("MAPQ difference (bhmem - yap)")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of MAPQ difference")

        # MAPQ distributions
        ax = axes[2]
        ax.hist(mapq_b, bins=43, alpha=0.5, label="bhmem", density=True)
        ax.hist(mapq_y, bins=43, alpha=0.5, label="yap", density=True)
        ax.set_xlabel("MAPQ")
        ax.set_ylabel("Density")
        ax.set_title("MAPQ distributions")
        ax.legend()

        plt.tight_layout()
        plot_path = f"{args.output_prefix}.plots.pdf"
        plt.savefig(plot_path)
        plt.close()
        print(f"Wrote plots to {plot_path}")
    elif args.plot and not HAS_MATPLOTLIB:
        print("WARNING: matplotlib not available, skipping plots.", file=sys.stderr)


if __name__ == "__main__":
    main()
