#!/usr/bin/env python
"""
Per-cell QC stats for single-cell HiC pipeline.
Designed to run as a SLURM array job (one task per cell).

Outputs: {output_dir}/{cell_id}.qc_stats.csv

Usage:
    python schic_qc_per_cell.py \
        --cell_id SRR1234567 \
        --fastp_dir trimmed_fastq \
        --bam_dir alignment \
        --output_dir qc_stats \
        --mapq 30
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
        "AdapterTrimmedReads": adapter.get("adapter_trimmed_reads", 0),
        "AdapterTrimmedRate": (
            adapter.get("adapter_trimmed_reads", 0) / input_reads * 100
            if input_reads > 0 else 0
        ),
        "Q20Rate_Before": before.get("q20_rate", 0) * 100,
        "Q30Rate_Before": before.get("q30_rate", 0) * 100,
        "Q20Rate_After": after.get("q20_rate", 0) * 100,
        "Q30Rate_After": after.get("q30_rate", 0) * 100,
        "GC_Before": before.get("gc_content", 0) * 100,
        "GC_After": after.get("gc_content", 0) * 100,
        "DuplicationRate_fastp": data.get("duplication", {}).get("rate", 0) * 100,
    }
    return stats


# ---------------------------------------------------------------------------
# 2. BAM stats via samtools flagstat (fast)
# ---------------------------------------------------------------------------

def parse_flagstat(bam_path):
    """Run samtools flagstat and parse output."""
    stats_text = pysam.flagstat(str(bam_path))
    results = {}
    for line in stats_text.strip().split("\n"):
        parts = line.split(" + ")
        primary_count = int(parts[0])
        rest = parts[1] if len(parts) > 1 else ""

        if "in total" in line:
            results["TotalAlignments"] = primary_count
        elif "secondary" in line:
            results["SecondaryAlignments"] = primary_count
        elif "supplementary" in line:
            results["SupplementaryAlignments"] = primary_count
        elif "mapped" in line and "mate mapped" not in line and "primary mapped" not in line:
            results["MappedAlignments"] = primary_count
        elif "properly paired" in line:
            results["ProperlyPaired"] = primary_count
        elif "duplicates" in line and "primary duplicates" not in line:
            results["Duplicates"] = primary_count
        elif "primary mapped" in line:
            results["PrimaryMapped"] = primary_count
        elif "primary duplicates" in line:
            results["PrimaryDuplicates"] = primary_count

    return results


def count_mapq_reads(bam_path, mapq_cutoff=30):
    """Count primary reads passing MAPQ cutoff using index stats for speed,
    falling back to iteration."""
    count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality >= mapq_cutoff:
                count += 1
    return count


# ---------------------------------------------------------------------------
# 3. Collect stats for one cell
# ---------------------------------------------------------------------------

def collect_cell_stats(cell_id, fastp_dir, bam_dir, mapq_cutoff=30):
    fastp_dir = pathlib.Path(fastp_dir)
    bam_dir = pathlib.Path(bam_dir)

    stats = {"CellID": cell_id}

    # --- fastp ---
    fastp_json = fastp_dir / f"{cell_id}.fastp.json"
    if fastp_json.exists():
        stats.update(parse_fastp_json(fastp_json))
    else:
        print(f"  [WARN] Missing fastp JSON: {fastp_json}")

    # --- sorted BAM (pre-dedup) ---
    sorted_bam = bam_dir / f"{cell_id}.sorted.bam"
    if sorted_bam.exists():
        flagstat = parse_flagstat(sorted_bam)
        stats["TotalAlignments"] = flagstat.get("TotalAlignments", 0)
        stats["MappedReads"] = flagstat.get("PrimaryMapped", flagstat.get("MappedAlignments", 0))
        stats["UnmappedReads"] = stats["TotalAlignments"] - stats["MappedReads"] - flagstat.get("SecondaryAlignments", 0) - flagstat.get("SupplementaryAlignments", 0)
        stats["SecondaryAlignments"] = flagstat.get("SecondaryAlignments", 0)
        stats["SupplementaryAlignments"] = flagstat.get("SupplementaryAlignments", 0)
        stats["ProperPairs"] = flagstat.get("ProperlyPaired", 0) // 2  # count pairs not reads

        print(f"  Counting MAPQ>={mapq_cutoff} reads in sorted BAM...")
        stats[f"MappedReadsMapQ{mapq_cutoff}"] = count_mapq_reads(sorted_bam, mapq_cutoff)

        # Mapping rates
        trimmed_total = stats.get("TrimmedReads_R1", 0) + stats.get("TrimmedReads_R2", 0)
        trimmed_pairs = stats.get("TrimmedReadPairs", 0)
        if trimmed_total > 0:
            stats["MappingRate"] = stats["MappedReads"] / trimmed_total * 100
            stats[f"MappingRateMapQ{mapq_cutoff}"] = (
                stats[f"MappedReadsMapQ{mapq_cutoff}"] / trimmed_total * 100
            )
        if trimmed_pairs > 0:
            stats["ProperPairRate"] = stats["ProperPairs"] / trimmed_pairs * 100
    else:
        print(f"  [WARN] Missing sorted BAM: {sorted_bam}")

    # --- calmd BAM (post-dedup) ---
    calmd_bam = bam_dir / f"{cell_id}.calmd.bam"
    if calmd_bam.exists():
        calmd_flagstat = parse_flagstat(calmd_bam)
        stats["DeduppedReads"] = calmd_flagstat.get(
            "PrimaryMapped", calmd_flagstat.get("MappedAlignments", 0)
        )

        print(f"  Counting MAPQ>={mapq_cutoff} reads in calmd BAM...")
        stats[f"DeduppedReadsMapQ{mapq_cutoff}"] = count_mapq_reads(calmd_bam, mapq_cutoff)

        pre_dedup = stats.get("MappedReads", 0)
        stats["DuplicatedReads"] = pre_dedup - stats["DeduppedReads"]
        stats["DuplicationRate"] = (
            stats["DuplicatedReads"] / pre_dedup * 100 if pre_dedup > 0 else 0
        )
    else:
        print(f"  [WARN] Missing calmd BAM: {calmd_bam}")

    # --- Derived ---
    if "DeduppedReads" in stats:
        stats["FinalReadPairs"] = stats["DeduppedReads"] // 2
    if "InputReadPairs" in stats and "FinalReadPairs" in stats:
        inp = stats["InputReadPairs"]
        stats["OverallPassRate"] = stats["FinalReadPairs"] / inp * 100 if inp > 0 else 0

    return stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Per-cell QC for scHiC pipeline (SLURM-friendly)"
    )
    parser.add_argument("--cell_id", required=True, help="Cell / sample ID")
    parser.add_argument("--fastp_dir", required=True, help="Dir with fastp JSON reports")
    parser.add_argument("--bam_dir", required=True, help="Dir with sorted and calmd BAMs")
    parser.add_argument("--output_dir", required=True, help="Dir for per-cell QC CSVs")
    parser.add_argument("--mapq", type=int, default=30, help="MAPQ cutoff (default: 30)")
    args = parser.parse_args()

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print(f"Processing cell: {args.cell_id}")
    stats = collect_cell_stats(
        args.cell_id,
        fastp_dir=args.fastp_dir,
        bam_dir=args.bam_dir,
        mapq_cutoff=args.mapq,
    )

    out_path = pathlib.Path(args.output_dir) / f"{args.cell_id}.qc_stats.csv"
    pd.Series(stats).to_frame().T.set_index("CellID").to_csv(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
