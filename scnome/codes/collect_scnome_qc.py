#!/usr/bin/env python
"""
Collect per-cell QC CSVs into a single scNOMe QC summary table.

Usage:
    python collect_scnome_qc.py \\
        --qc_dir qc_stats \\
        --acc_list acc_list.txt \\
        --output scnome_qc_summary.csv
"""

import argparse
import pathlib
import sys

import pandas as pd


PRIORITY_COLS = [
    # Combined summary
    "TotalInputReads", "Combined_TotalReads", "Combined_UniqMapped",
    "Combined_MappingRate",
    "Combined_mCG_Rate", "Combined_mCHG_Rate", "Combined_mCHH_Rate",
    # Trim stats — R1
    "Trim_R1_InputReads", "Trim_R1_TrimmedReads", "Trim_R1_TrimmedReads_Pct",
    "Trim_R1_ReadsWithAdapters", "Trim_R1_ReadsWithAdapters_Pct",
    # Trim stats — R2
    "Trim_R2_InputReads", "Trim_R2_TrimmedReads", "Trim_R2_TrimmedReads_Pct",
    "Trim_R2_ReadsWithAdapters", "Trim_R2_ReadsWithAdapters_Pct",
    # Bismark — R1
    "Bismark_R1_TotalReads", "Bismark_R1_UniqMapped", "Bismark_R1_MappingRate",
    "Bismark_R1_Unmapped", "Bismark_R1_Ambiguous",
    "Bismark_R1_mCG_Rate", "Bismark_R1_mCHG_Rate", "Bismark_R1_mCHH_Rate",
    "Bismark_R1_TotalC",
    # Bismark — R2
    "Bismark_R2_TotalReads", "Bismark_R2_UniqMapped", "Bismark_R2_MappingRate",
    "Bismark_R2_Unmapped", "Bismark_R2_Ambiguous",
    "Bismark_R2_mCG_Rate", "Bismark_R2_mCHG_Rate", "Bismark_R2_mCHH_Rate",
    "Bismark_R2_TotalC",
    # BAM summary — R1
    "BAM_R1_TotalReads", "BAM_R1_PrimaryMappedReads",
    "BAM_R1_UniqMappedMapQ30Reads", "BAM_R1_DuplicatePrimaryMapped",
    "BAM_R1_UniqMappedMapQ30ReadsToLambda",
    "BAM_R1_UniqMappedMapQ30ReadsToTargetSpecies",
    # BAM summary — R2
    "BAM_R2_TotalReads", "BAM_R2_PrimaryMappedReads",
    "BAM_R2_UniqMappedMapQ30Reads", "BAM_R2_DuplicatePrimaryMapped",
    "BAM_R2_UniqMappedMapQ30ReadsToLambda",
    "BAM_R2_UniqMappedMapQ30ReadsToTargetSpecies",
]

OVERVIEW_COLS = [
    "TotalInputReads",
    "Combined_MappingRate",
    "Combined_mCG_Rate",
    "Combined_mCHG_Rate",
    "Combined_mCHH_Rate",
    "Bismark_R1_MappingRate",
    "Bismark_R2_MappingRate",
    "BAM_R1_DuplicatePrimaryMapped_Pct",
    "BAM_R2_DuplicatePrimaryMapped_Pct",
]


def main():
    parser = argparse.ArgumentParser(
        description="Collect per-cell scNOMe QC stats into one summary CSV"
    )
    parser.add_argument("--qc_dir", required=True,
                        help="Directory containing per-cell .qc_stats.csv files")
    parser.add_argument("--acc_list", required=True,
                        help="Text file with one cell ID per line")
    parser.add_argument("--output", default="scnome_qc_summary.csv",
                        help="Output CSV path")
    args = parser.parse_args()

    qc_dir = pathlib.Path(args.qc_dir)

    with open(args.acc_list) as f:
        expected = [line.strip() for line in f if line.strip()]

    found = {}
    for cell_id in expected:
        p = qc_dir / f"{cell_id}.qc_stats.csv"
        if p.exists():
            found[cell_id] = p

    missing = [c for c in expected if c not in found]
    if missing:
        print(f"WARNING: {len(missing)} cells missing QC stats:")
        for m in missing:
            print(f"  - {m}")

    if not found:
        print("ERROR: No .qc_stats.csv files found.")
        sys.exit(1)

    print(f"Found {len(found)} / {len(expected)} per-cell QC files")

    dfs = []
    for cell_id, path in found.items():
        try:
            dfs.append(pd.read_csv(path, index_col=0))
        except Exception as e:
            print(f"  [WARN] Failed to read {path}: {e}")

    if not dfs:
        print("ERROR: No valid QC files could be read.")
        sys.exit(1)

    summary = pd.concat(dfs).sort_index()
    summary.index.name = "CellID"

    existing = [c for c in PRIORITY_COLS if c in summary.columns]
    extra = [c for c in summary.columns if c not in PRIORITY_COLS]
    summary = summary[existing + extra]

    summary.to_csv(args.output)
    print(f"\nSaved: {args.output}")
    print(f"Shape: {summary.shape[0]} cells x {summary.shape[1]} columns")

    print("\n--- Key Metrics Overview ---")
    for col in OVERVIEW_COLS:
        if col in summary.columns:
            s = pd.to_numeric(summary[col], errors="coerce").dropna()
            if len(s) > 0:
                print(f"  {col:50s}  median={s.median():>9.2f}  "
                      f"mean={s.mean():>9.2f}  "
                      f"range=[{s.min():.2f}, {s.max():.2f}]")

    print("\nDone.")


if __name__ == "__main__":
    main()
