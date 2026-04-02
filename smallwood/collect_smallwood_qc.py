#!/usr/bin/env python
"""
Collect per-cell QC CSVs into a single summary table.

Usage:
    python collect_smallwood_qc.py \
        --qc_dir qc_stats \
        --acc_list acc_list.txt \
        --output smallwood_qc_summary.csv
"""

import argparse
import pathlib
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Collect per-cell QC stats into one summary CSV"
    )
    parser.add_argument("--qc_dir", required=True)
    parser.add_argument("--acc_list", required=True)
    parser.add_argument("--output", default="smallwood_qc_summary.csv")
    args = parser.parse_args()

    qc_dir = pathlib.Path(args.qc_dir)

    with open(args.acc_list) as f:
        expected = [line.strip() for line in f if line.strip()]

    # Find all per-cell CSVs
    found_files = {}
    for cell_id in expected:
        path = qc_dir / f"{cell_id}.qc_stats.csv"
        if path.exists():
            found_files[cell_id] = path

    missing = [c for c in expected if c not in found_files]
    if missing:
        print(f"WARNING: {len(missing)} cells missing QC stats:")
        for m in missing:
            print(f"  - {m}")

    if not found_files:
        print("ERROR: No .qc_stats.csv files found.")
        sys.exit(1)

    print(f"Found {len(found_files)} / {len(expected)} per-cell QC files")

    # Read and concatenate
    dfs = []
    for cell_id, path in found_files.items():
        try:
            df = pd.read_csv(path, index_col=0)
            dfs.append(df)
        except Exception as e:
            print(f"  [WARN] Failed to read {path}: {e}")

    summary = pd.concat(dfs)
    summary.index.name = "CellID"
    summary = summary.sort_index()

    # Reorder columns
    priority_cols = [
        "Trim_R1_InputReads", "Trim_R1_TrimmedReads", "Trim_R1_TrimmedReads_Pct",
        "Trim_R2_InputReads", "Trim_R2_TrimmedReads", "Trim_R2_TrimmedReads_Pct",
        "Bismark_TotalReads", "Bismark_UniqMapped", "Bismark_MappingRate",
        "Bismark_Unmapped", "Bismark_Ambiguous",
        "Bismark_mCG_Rate", "Bismark_mCHG_Rate", "Bismark_mCHH_Rate",
        "Bismark_TotalC",
        "BAM_TotalAlignments", "BAM_PrimaryMapped", "BAM_Duplicates",
        "BAM_MapQ30", "BAM_MapQ30_Rate",
        "chrM_noncpg", "chrM_endo", "chrM_exo",
        "chr19_noncpg", "chr19_endo", "chr19_exo",
        "CpG_TotalSites", "CpG_TotalReads", "CpG_MeanMethylation",
        "CpG_MeanCoverage", "CpG_MethylatedSites", "CpG_MethylatedFrac",
    ]
    existing = [c for c in priority_cols if c in summary.columns]
    extra = [c for c in summary.columns if c not in priority_cols]
    summary = summary[existing + extra]

    summary.to_csv(args.output)
    print(f"\nSaved: {args.output}")
    print(f"Shape: {summary.shape[0]} cells x {summary.shape[1]} columns")

    # Print overview
    print("\n--- Key Metrics Overview ---")
    overview_cols = [
        "Trim_R1_InputReads", "Trim_R1_TrimmedReads_Pct",
        "Bismark_MappingRate", "Bismark_mCG_Rate",
        "BAM_MapQ30", "BAM_MapQ30_Rate",
        "chrM_noncpg", "chrM_endo", "chrM_exo",
        "chr19_noncpg", "chr19_endo", "chr19_exo",
        "CpG_TotalSites", "CpG_TotalReads", "CpG_MeanMethylation",
    ]
    for col in overview_cols:
        if col in summary.columns:
            s = pd.to_numeric(summary[col], errors="coerce").dropna()
            if len(s) > 0:
                print(f"  {col:40s}  median={s.median():>12.2f}  "
                      f"mean={s.mean():>12.2f}  "
                      f"range=[{s.min():.2f}, {s.max():.2f}]")

    print("\nDone.")


if __name__ == "__main__":
    main()
