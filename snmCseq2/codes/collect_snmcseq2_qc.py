#!/usr/bin/env python
"""
Collect per-cell QC CSVs into a single snmC-seq2 QC summary table.

Usage:
    python collect_snmcseq2_qc.py \\
        --qc_dir qc_stats \\
        --acc_list acc_list.txt \\
        --output snmcseq2_qc_summary.csv
"""

import argparse
import pathlib
import sys

import pandas as pd


PRIORITY_COLS = [
    # Metadata
    "Genome",
    # Combined summary
    "Combined_TotalReads", "Combined_UniqMapped", "Combined_MappingRate",
    "Combined_mCG_Rate", "Combined_mCHG_Rate", "Combined_mCHH_Rate",
    "BAM_Combined_Mapped", "BAM_Combined_MappedMapQ30",
    # Bismark — R1 (PBAT)
    "Bismark_R1_TotalReads", "Bismark_R1_UniqMapped", "Bismark_R1_MappingRate",
    "Bismark_R1_Unmapped", "Bismark_R1_Ambiguous",
    "Bismark_R1_mCG_Rate", "Bismark_R1_mCHG_Rate", "Bismark_R1_mCHH_Rate",
    "Bismark_R1_TotalC",
    # Bismark — R2 (non-directional)
    "Bismark_R2_TotalReads", "Bismark_R2_UniqMapped", "Bismark_R2_MappingRate",
    "Bismark_R2_Unmapped", "Bismark_R2_Ambiguous",
    "Bismark_R2_mCG_Rate", "Bismark_R2_mCHG_Rate", "Bismark_R2_mCHH_Rate",
    "Bismark_R2_TotalC",
    # BAM summary — R1
    "BAM_R1_TotalReads", "BAM_R1_Mapped", "BAM_R1_MappedMapQ30",
    "BAM_R1_MappedMapQ30ToLambda", "BAM_R1_MappedMapQ30ToTargetSpecies",
    "BAM_R1_DuplicateMapped",
    # BAM summary — R2
    "BAM_R2_TotalReads", "BAM_R2_Mapped", "BAM_R2_MappedMapQ30",
    "BAM_R2_MappedMapQ30ToLambda", "BAM_R2_MappedMapQ30ToTargetSpecies",
    "BAM_R2_DuplicateMapped",
]

OVERVIEW_COLS = [
    "Combined_TotalReads",
    "Combined_MappingRate",
    "Combined_mCG_Rate",
    "Combined_mCHG_Rate",
    "Combined_mCHH_Rate",
    "Bismark_R1_MappingRate",
    "Bismark_R2_MappingRate",
    "BAM_R1_MappedMapQ30_Pct",
    "BAM_R2_MappedMapQ30_Pct",
    "BAM_R1_MappedMapQ30ToTargetSpecies",
    "BAM_R2_MappedMapQ30ToTargetSpecies",
]


def main():
    parser = argparse.ArgumentParser(
        description="Collect per-cell snmC-seq2 QC stats into one summary CSV"
    )
    parser.add_argument("--qc_dir", required=True,
                        help="Directory containing per-cell .qc_stats.csv files")
    parser.add_argument("--acc_list", required=True,
                        help="Text file with one cell ID per line")
    parser.add_argument("--output", default="snmcseq2_qc_summary.csv",
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

    # Per-genome breakdown if Genome column present
    if "Genome" in summary.columns:
        for genome, grp in summary.groupby("Genome"):
            print(f"\n  Genome={genome}: {len(grp)} cells")

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
