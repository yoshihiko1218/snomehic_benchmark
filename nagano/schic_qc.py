#!/usr/bin/env python
"""
QC Report Generator for Single-Cell HiC Pipeline
=================================================

Collects QC stats from:
  1. fastp trimming (JSON reports)
  2. Bowtie2 alignment (sorted BAM)
  3. Post-dedup BAM (calmd BAM)

Produces a per-cell QC summary CSV similar to the m3c MappingSummary.

Usage:
    python schic_qc.py \
        --acc_list acc_list.txt \
        --fastp_dir trimmed_fastq \
        --bam_dir alignment \
        --output qc_summary.csv

Dependencies:
    - pysam
    - pandas
    - json (stdlib)
"""

import argparse
import json
import pathlib
from collections import defaultdict

import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# 1. Parse fastp JSON
# ---------------------------------------------------------------------------

def parse_fastp_json(json_path):
    """
    Extract key QC metrics from a fastp JSON report.

    Returns a dict with:
        - InputReadPairs, InputReads_R1, InputReads_R2
        - InputBP_R1, InputBP_R2
        - TrimmedReadPairs, TrimmedReads_R1, TrimmedReads_R2
        - TrimmedBP_R1, TrimmedBP_R2
        - TrimmedReadsRate
        - AdapterTrimmedReads
        - Q20Rate_Before, Q30Rate_Before
        - Q20Rate_After, Q30Rate_After
        - GC_Before, GC_After
        - DuplicationRate_fastp (fastp's own estimate)
    """
    with open(json_path) as f:
        data = json.load(f)

    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    adapter = data.get("adapter_cutting", {})
    read1_before = data.get("read1_before_filtering", {})
    read2_before = data.get("read2_before_filtering", {})
    read1_after = data.get("read1_after_filtering", {})
    read2_after = data.get("read2_after_filtering", {})

    input_reads = before.get("total_reads", 0)
    output_reads = after.get("total_reads", 0)

    stats = {
        # Read counts
        "InputReadPairs": input_reads // 2,
        "InputReads_R1": read1_before.get("total_reads", input_reads // 2),
        "InputReads_R2": read2_before.get("total_reads", input_reads // 2),
        "InputBP_R1": read1_before.get("total_bases", 0),
        "InputBP_R2": read2_before.get("total_bases", 0),
        "TrimmedReadPairs": output_reads // 2,
        "TrimmedReads_R1": read1_after.get("total_reads", output_reads // 2),
        "TrimmedReads_R2": read2_after.get("total_reads", output_reads // 2),
        "TrimmedBP_R1": read1_after.get("total_bases", 0),
        "TrimmedBP_R2": read2_after.get("total_bases", 0),
        "TrimmedReadsRate": output_reads / input_reads * 100 if input_reads > 0 else 0,

        # Adapter
        "AdapterTrimmedReads": adapter.get("adapter_trimmed_reads", 0),
        "AdapterTrimmedRate": (
            adapter.get("adapter_trimmed_reads", 0) / input_reads * 100
            if input_reads > 0 else 0
        ),

        # Quality
        "Q20Rate_Before": before.get("q20_rate", 0) * 100,
        "Q30Rate_Before": before.get("q30_rate", 0) * 100,
        "Q20Rate_After": after.get("q20_rate", 0) * 100,
        "Q30Rate_After": after.get("q30_rate", 0) * 100,

        # GC
        "GC_Before": before.get("gc_content", 0) * 100,
        "GC_After": after.get("gc_content", 0) * 100,

        # fastp duplication estimate
        "DuplicationRate_fastp": data.get("duplication", {}).get("rate", 0) * 100,
    }
    return stats


# ---------------------------------------------------------------------------
# 2. BAM-level stats
# ---------------------------------------------------------------------------

def count_bam_stats(bam_path, mapq_cutoff=30):
    """
    Count alignment-level stats from a coordinate-sorted, indexed BAM.

    Returns a dict with:
        - TotalAlignments
        - MappedReads (not unmapped, not secondary, not supplementary)
        - UnmappedReads
        - MappedReadsMapQ{mapq_cutoff}
        - ProperPairs (both mates properly paired, counted once per pair)
        - SecondaryAlignments
        - SupplementaryAlignments
    """
    stats = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            stats["TotalAlignments"] += 1

            if read.is_unmapped:
                stats["UnmappedReads"] += 1
                continue

            if read.is_secondary:
                stats["SecondaryAlignments"] += 1
                continue
            if read.is_supplementary:
                stats["SupplementaryAlignments"] += 1
                continue

            # Primary alignment
            stats["MappedReads"] += 1

            if read.mapping_quality >= mapq_cutoff:
                stats[f"MappedReadsMapQ{mapq_cutoff}"] += 1

            if read.is_proper_pair and read.is_read1:
                stats["ProperPairs"] += 1

    return dict(stats)


def count_dedup_stats(sorted_bam_path, calmd_bam_path, mapq_cutoff=30):
    """
    Compare pre-dedup (sorted) and post-dedup (calmd) BAMs to derive
    deduplication statistics.

    Returns a dict with:
        - DeduppedReads (primary mapped in calmd BAM)
        - DeduppedReadsMapQ{mapq_cutoff}
        - DuplicatedReads
        - DuplicationRate (%)
    """
    def _count_primary(bam_path, mapq_cutoff):
        total = 0
        mapq_pass = 0
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                total += 1
                if read.mapping_quality >= mapq_cutoff:
                    mapq_pass += 1
        return total, mapq_pass

    pre_total, pre_mapq = _count_primary(sorted_bam_path, mapq_cutoff)
    post_total, post_mapq = _count_primary(calmd_bam_path, mapq_cutoff)

    duplicated = pre_total - post_total
    dup_rate = duplicated / pre_total * 100 if pre_total > 0 else 0

    return {
        "PreDedupReads": pre_total,
        "DeduppedReads": post_total,
        f"DeduppedReadsMapQ{mapq_cutoff}": post_mapq,
        "DuplicatedReads": duplicated,
        "DuplicationRate": dup_rate,
    }


# ---------------------------------------------------------------------------
# 3. Per-cell summary
# ---------------------------------------------------------------------------

def collect_cell_stats(
    cell_id,
    fastp_dir,
    bam_dir,
    mapq_cutoff=30,
):
    """Collect all QC stats for one cell."""
    fastp_dir = pathlib.Path(fastp_dir)
    bam_dir = pathlib.Path(bam_dir)

    stats = {"CellID": cell_id}

    # --- fastp ---
    fastp_json = fastp_dir / f"{cell_id}.fastp.json"
    if fastp_json.exists():
        stats.update(parse_fastp_json(fastp_json))
    else:
        print(f"  [WARN] Missing fastp JSON: {fastp_json}")

    # --- sorted BAM (pre-dedup alignment stats) ---
    sorted_bam = bam_dir / f"{cell_id}.sorted.bam"
    if sorted_bam.exists():
        bam_stats = count_bam_stats(sorted_bam, mapq_cutoff=mapq_cutoff)
        stats.update(bam_stats)

        # Mapping rate
        trimmed_pairs = stats.get("TrimmedReadPairs", 0)
        total_trimmed_reads = trimmed_pairs * 2 if trimmed_pairs else 0
        if total_trimmed_reads > 0:
            stats["MappingRate"] = (
                bam_stats.get("MappedReads", 0) / total_trimmed_reads * 100
            )
            stats[f"MappingRateMapQ{mapq_cutoff}"] = (
                bam_stats.get(f"MappedReadsMapQ{mapq_cutoff}", 0)
                / total_trimmed_reads * 100
            )
            stats["ProperPairRate"] = (
                bam_stats.get("ProperPairs", 0) / trimmed_pairs * 100
                if trimmed_pairs > 0 else 0
            )
    else:
        print(f"  [WARN] Missing sorted BAM: {sorted_bam}")

    # --- calmd BAM (post-dedup) ---
    calmd_bam = bam_dir / f"{cell_id}.calmd.bam"
    if sorted_bam.exists() and calmd_bam.exists():
        dedup_stats = count_dedup_stats(
            sorted_bam, calmd_bam, mapq_cutoff=mapq_cutoff
        )
        stats.update(dedup_stats)
    else:
        if not calmd_bam.exists():
            print(f"  [WARN] Missing calmd BAM: {calmd_bam}")

    return stats


# ---------------------------------------------------------------------------
# 4. Aggregate + derived columns
# ---------------------------------------------------------------------------

def add_derived_columns(df):
    """Add summary / derived QC columns to the dataframe."""
    # Final usable reads (post-dedup)
    if "DeduppedReads" in df.columns:
        df["FinalReadPairs"] = df["DeduppedReads"] // 2

    # Overall pass rate: from raw input to final deduped
    if "InputReadPairs" in df.columns and "FinalReadPairs" in df.columns:
        df["OverallPassRate"] = (
            df["FinalReadPairs"] / df["InputReadPairs"] * 100
        ).where(df["InputReadPairs"] > 0, 0)

    return df


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate QC summary for single-cell HiC pipeline"
    )
    parser.add_argument(
        "--acc_list", required=True,
        help="Text file with one cell/sample ID per line",
    )
    parser.add_argument(
        "--fastp_dir", required=True,
        help="Directory containing fastp JSON reports",
    )
    parser.add_argument(
        "--bam_dir", required=True,
        help="Directory containing sorted and calmd BAM files",
    )
    parser.add_argument(
        "--output", default="qc_summary.csv",
        help="Output CSV path (default: qc_summary.csv)",
    )
    parser.add_argument(
        "--mapq", type=int, default=30,
        help="MAPQ cutoff for high-quality reads (default: 30)",
    )
    args = parser.parse_args()

    # Read cell IDs
    with open(args.acc_list) as f:
        cell_ids = [line.strip() for line in f if line.strip()]

    print(f"Collecting QC stats for {len(cell_ids)} cells...")

    records = []
    for cell_id in cell_ids:
        print(f"  Processing {cell_id}")
        stats = collect_cell_stats(
            cell_id,
            fastp_dir=args.fastp_dir,
            bam_dir=args.bam_dir,
            mapq_cutoff=args.mapq,
        )
        records.append(stats)

    df = pd.DataFrame(records).set_index("CellID")
    df = add_derived_columns(df)

    # Sort columns into logical groups
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
        f"MappedReadsMapQ{args.mapq}",
        "SecondaryAlignments", "SupplementaryAlignments",
        "ProperPairs",
        "MappingRate", f"MappingRateMapQ{args.mapq}", "ProperPairRate",
        # Dedup
        "PreDedupReads", "DeduppedReads", f"DeduppedReadsMapQ{args.mapq}",
        "DuplicatedReads", "DuplicationRate",
        # Summary
        "FinalReadPairs", "OverallPassRate",
    ]
    existing = [c for c in col_order if c in df.columns]
    extra = [c for c in df.columns if c not in col_order]
    df = df[existing + extra]

    df.to_csv(args.output)
    print(f"\nQC summary saved to {args.output}")
    print(f"Shape: {df.shape[0]} cells x {df.shape[1]} metrics")

    # Print a quick overview
    print("\n--- Quick Overview ---")
    overview_cols = [
        "InputReadPairs", "TrimmedReadsRate", "MappingRate",
        "DuplicationRate", "FinalReadPairs", "OverallPassRate",
    ]
    for col in overview_cols:
        if col in df.columns:
            print(f"  {col}: median={df[col].median():.1f}, "
                  f"mean={df[col].mean():.1f}, "
                  f"min={df[col].min():.1f}, max={df[col].max():.1f}")


if __name__ == "__main__":
    main()
