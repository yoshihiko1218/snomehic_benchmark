#!/usr/bin/env python
"""
Collect per-cell QC CSVs into a single summary table.

Usage:
    python collect_qc.py \
        --qc_dir qc_stats \
        --output qc_summary.csv \
        [--acc_list acc_list.txt]

If --acc_list is provided, warns about any missing cells.
"""

import argparse
import pathlib
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Collect per-cell QC stats into one summary CSV"
    )
    parser.add_argument(
        "--qc_dir", required=True,
        help="Directory containing per-cell .qc_stats.csv files",
    )
    parser.add_argument(
        "--output", default="qc_summary.csv",
        help="Output combined CSV (default: qc_summary.csv)",
    )
    parser.add_argument(
        "--acc_list", default=None,
        help="Optional: acc_list.txt to check for missing cells",
    )
    args = parser.parse_args()

    qc_dir = pathlib.Path(args.qc_dir)
    csv_files = sorted(qc_dir.glob("*.qc_stats.csv"))

    if len(csv_files) == 0:
        print(f"ERROR: No .qc_stats.csv files found in {qc_dir}")
        sys.exit(1)

    print(f"Found {len(csv_files)} per-cell QC files")

    # Check for missing cells
    if args.acc_list:
        with open(args.acc_list) as f:
            expected = {line.strip() for line in f if line.strip()}
        found = {p.name.replace(".qc_stats.csv", "") for p in csv_files}
        missing = expected - found
        if missing:
            print(f"\nWARNING: {len(missing)} cells missing QC stats:")
            for m in sorted(missing):
                print(f"  - {m}")
            print()

    # Read and concatenate
    dfs = []
    for path in csv_files:
        try:
            df = pd.read_csv(path, index_col=0)
            dfs.append(df)
        except Exception as e:
            print(f"  [WARN] Failed to read {path}: {e}")

    if not dfs:
        print("ERROR: No valid QC files could be read")
        sys.exit(1)

    summary = pd.concat(dfs)
    summary.index.name = "CellID"
    summary = summary.sort_index()

    # Reorder columns into logical groups
    col_order = [
        # Trimming
        "InputReadPairs", "InputReads_R1", "InputReads_R2",
        "InputBP_R1", "InputBP_R2",
        "TrimmedReadPairs", "TrimmedReads_R1", "TrimmedReads_R2",
        "TrimmedBP_R1", "TrimmedBP_R2",
        "TrimmedReadsRate",
        "AdapterTrimmedReads", "AdapterTrimmedRate",
        # Quality
        "Q20Rate_Before", "Q30Rate_Before",
        "Q20Rate_After", "Q30Rate_After",
        "GC_Before", "GC_After",
        "DuplicationRate_fastp",
        # Alignment
        "TotalAlignments", "MappedReads", "UnmappedReads",
        "MappedReadsMapQ30",
        "SecondaryAlignments", "SupplementaryAlignments",
        "ProperPairs",
        "MappingRate", "MappingRateMapQ30", "ProperPairRate",
        # Dedup
        "DeduppedReads", "DeduppedReadsMapQ30",
        "DuplicatedReads", "DuplicationRate",
        # Summary
        "FinalReadPairs", "OverallPassRate",
    ]
    existing = [c for c in col_order if c in summary.columns]
    extra = [c for c in summary.columns if c not in col_order]
    summary = summary[existing + extra]

    summary.to_csv(args.output)
    print(f"\nSaved combined QC summary: {args.output}")
    print(f"Shape: {summary.shape[0]} cells x {summary.shape[1]} metrics")

    # Quick overview
    print("\n--- Quick Overview ---")
    overview_cols = [
        "InputReadPairs", "TrimmedReadsRate", "MappingRate",
        "MappingRateMapQ30", "ProperPairRate",
        "DuplicationRate", "FinalReadPairs", "OverallPassRate",
    ]
    for col in overview_cols:
        if col in summary.columns:
            s = summary[col].dropna()
            if len(s) > 0:
                print(f"  {col:30s}  median={s.median():>12.1f}  "
                      f"mean={s.mean():>12.1f}  "
                      f"min={s.min():>12.1f}  max={s.max():>12.1f}")


if __name__ == "__main__":
    main()
