#!/usr/bin/env python
"""
Per-cell QC for Rupture Droplet Hi-C Pipeline
==============================================

Generates per-barcode (single-cell) QC stats directly from pipeline outputs.
No dependency on any pre-existing QC files.

Reads from:
  - 03.mapping/{sample}_{genome}.bam             (CB-tagged BAM, name-sorted stream)
  - 03.mapping/{sample}_{genome}.sc.pairs.gz     (deduplicated pairs with CB column)

Produces:
  - {output_prefix}.per_cell.tsv       — one row per barcode, all metrics
  - {output_prefix}.sample_summary.tsv — one row per sample, aggregate stats

Usage:
    python rupture_qc.py \\
        --bam 03.mapping/SRR27586278_hg38.bam \\
        --pairs 03.mapping/SRR27586278_hg38.sc.pairs.gz \\
        --output SRR27586278_hg38_qc \\
        --mapq 30

    # Or just BAM (if pairs not yet available):
    python rupture_qc.py \\
        --bam 03.mapping/SRR27586278_hg38.bam \\
        --output SRR27586278_hg38_qc

    # Or just pairs:
    python rupture_qc.py \\
        --pairs 03.mapping/SRR27586278_hg38.sc.pairs.gz \\
        --output SRR27586278_hg38_qc

Dependencies: pysam, pandas
"""

import argparse
import gzip
import pathlib
from collections import defaultdict

import pandas as pd
import pysam


# ═══════════════════════════════════════════════════════════════════
# 1. BAM-based per-cell stats
# ═══════════════════════════════════════════════════════════════════

def get_cb(read):
    """Extract cell barcode from CB tag or read name."""
    if read.has_tag("CB"):
        return read.get_tag("CB")
    qn = read.query_name
    return qn.split(":")[-1] if ":" in qn else "NO_CB"


def qc_from_bam(bam_path, mapq_cutoff=30):
    """
    Iterate through a BAM (name-sorted stream via pysam) and collect
    per-barcode fragment-level stats.

    Pairs consecutive reads by query name. For each fragment (read pair):
      - TotalFragments: all fragments
      - UniqMapped: both mates primary + mapped
      - UniqMappedMapQ{mapq}: both mates pass MAPQ
      - Duplicates: either mate flagged as duplicate
      - UniqMappedMapQNoDup: mapped, pass MAPQ, not duplicate
      - UniqMappedMapQNoDup_Cis / Trans: same/different chromosome
      - Cis_1kb, Cis_10kb, Cis_20kb: cis contacts at distance thresholds

    Returns dict[barcode] -> dict[metric] -> int
    """
    print(f"  Reading BAM: {bam_path}")
    counts = defaultdict(lambda: defaultdict(int))

    mq_label = f"MapQ{mapq_cutoff}"

    mode = "rb" if str(bam_path).endswith(".bam") else "rc"
    with pysam.AlignmentFile(str(bam_path), mode) as bam:
        prev_read = None
        n_reads = 0

        for read in bam.fetch(until_eof=True):
            n_reads += 1
            if n_reads % 10_000_000 == 0:
                print(f"    processed {n_reads:,} alignments, {len(counts):,} barcodes...")

            if prev_read is None:
                prev_read = read
                continue

            if prev_read.query_name != read.query_name:
                prev_read = read
                continue

            r1, r2 = prev_read, read
            prev_read = None

            cb = get_cb(r1)
            c = counts[cb]
            c["TotalFragments"] += 1

            # Both mates primary and mapped
            if (r1.is_paired
                    and not r1.is_unmapped and not r2.is_unmapped
                    and not r1.is_secondary and not r2.is_secondary
                    and not r1.is_supplementary and not r2.is_supplementary):
                c["UniqMapped"] += 1

                if r1.mapping_quality >= mapq_cutoff and r2.mapping_quality >= mapq_cutoff:
                    c[f"UniqMapped{mq_label}"] += 1

                    is_dup = r1.is_duplicate or r2.is_duplicate
                    if is_dup:
                        c["Duplicates"] += 1
                    else:
                        c[f"UniqMapped{mq_label}NoDup"] += 1

                        if r1.reference_name != r1.next_reference_name:
                            c[f"UniqMapped{mq_label}NoDup_Trans"] += 1
                        else:
                            c[f"UniqMapped{mq_label}NoDup_Cis"] += 1
                            tlen = abs(r1.template_length)
                            if tlen >= 1000:
                                c[f"UniqMapped{mq_label}NoDup_Cis1kb"] += 1
                            if tlen >= 10000:
                                c[f"UniqMapped{mq_label}NoDup_Cis10kb"] += 1
                            if tlen >= 20000:
                                c[f"UniqMapped{mq_label}NoDup_Cis20kb"] += 1

    print(f"    done. {n_reads:,} alignments, {len(counts):,} barcodes")
    return dict(counts)


# ═══════════════════════════════════════════════════════════════════
# 2. Pairs-based per-cell stats
# ═══════════════════════════════════════════════════════════════════

def find_cb_column(header_fields):
    """Find the CB column index in a pairs file header."""
    for i, field in enumerate(header_fields):
        if field.upper() in ("CB", "CB1"):
            return i
    return None


def qc_from_pairs(pairs_path):
    """
    Parse a .pairs.gz file (pairtools output with CB column) and collect
    per-barcode contact stats.

    For each pair:
      - TotalPairs
      - Cis / Trans
      - Cis_1kb, Cis_10kb, Cis_20kb (by genomic distance)
      - pair_type counts (UU, UR, RU, etc.)

    Returns dict[barcode] -> dict[metric] -> int
    """
    print(f"  Reading pairs: {pairs_path}")
    counts = defaultdict(lambda: defaultdict(int))

    opener = gzip.open if str(pairs_path).endswith(".gz") else open

    cb_col = None
    chrom1_col = 1  # standard pairs columns
    pos1_col = 2
    chrom2_col = 3
    pos2_col = 4
    pair_type_col = 7

    n_pairs = 0

    with opener(str(pairs_path), "rt") as f:
        for line in f:
            if line.startswith("#columns:"):
                fields = line.strip().split(":")[1].strip().split()
                # Map column names to indices
                col_map = {name: i for i, name in enumerate(fields)}
                chrom1_col = col_map.get("chrom1", 1)
                pos1_col = col_map.get("pos1", 2)
                chrom2_col = col_map.get("chrom2", 3)
                pos2_col = col_map.get("pos2", 4)
                pair_type_col = col_map.get("pair_type", 7)
                cb_col = col_map.get("CB", col_map.get("CB1", None))
                continue
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            n_pairs += 1

            if n_pairs % 10_000_000 == 0:
                print(f"    processed {n_pairs:,} pairs, {len(counts):,} barcodes...")

            # Get barcode
            if cb_col is not None and cb_col < len(parts):
                cb = parts[cb_col]
            else:
                cb = "NO_CB"

            c = counts[cb]
            c["TotalPairs"] += 1

            chrom1 = parts[chrom1_col]
            chrom2 = parts[chrom2_col]

            # Pair type
            if pair_type_col < len(parts):
                pt = parts[pair_type_col]
                c[f"PairType_{pt}"] += 1

            # Cis / Trans
            if chrom1 == chrom2:
                c["Pairs_Cis"] += 1
                try:
                    dist = abs(int(parts[pos2_col]) - int(parts[pos1_col]))
                    if dist >= 1000:
                        c["Pairs_Cis1kb"] += 1
                    if dist >= 10000:
                        c["Pairs_Cis10kb"] += 1
                    if dist >= 20000:
                        c["Pairs_Cis20kb"] += 1
                except (ValueError, IndexError):
                    pass
            else:
                c["Pairs_Trans"] += 1

    print(f"    done. {n_pairs:,} pairs, {len(counts):,} barcodes")
    return dict(counts)


# ═══════════════════════════════════════════════════════════════════
# 3. Merge and compute derived stats
# ═══════════════════════════════════════════════════════════════════

def merge_stats(bam_counts, pairs_counts):
    """Merge BAM-based and pairs-based per-cell dicts."""
    all_barcodes = set()
    if bam_counts:
        all_barcodes.update(bam_counts.keys())
    if pairs_counts:
        all_barcodes.update(pairs_counts.keys())

    merged = {}
    for cb in all_barcodes:
        d = {}
        if bam_counts and cb in bam_counts:
            d.update(bam_counts[cb])
        if pairs_counts and cb in pairs_counts:
            d.update(pairs_counts[cb])
        merged[cb] = d
    return merged


def build_dataframe(merged_counts, mapq_cutoff=30):
    """Convert merged per-cell dict to DataFrame with derived columns."""
    df = pd.DataFrame.from_dict(merged_counts, orient="index")
    df.index.name = "CB"
    df = df.fillna(0).astype(int, errors="ignore")
    df = df.sort_index()

    mq = f"MapQ{mapq_cutoff}"

    # ── Derived ratios (BAM-based) ──
    if "TotalFragments" in df.columns and "UniqMapped" in df.columns:
        df["MappingRate"] = (
            df["UniqMapped"] / df["TotalFragments"]
        ).where(df["TotalFragments"] > 0, 0)

    mq_col = f"UniqMapped{mq}"
    if mq_col in df.columns and "TotalFragments" in df.columns:
        df[f"{mq}_Rate"] = (
            df[mq_col] / df["TotalFragments"]
        ).where(df["TotalFragments"] > 0, 0)

    nodup_col = f"UniqMapped{mq}NoDup"
    if nodup_col in df.columns and mq_col in df.columns:
        df["DuplicationRate"] = (
            1 - df[nodup_col] / df[mq_col]
        ).where(df[mq_col] > 0, 0)

    cis_col = f"UniqMapped{mq}NoDup_Cis"
    trans_col = f"UniqMapped{mq}NoDup_Trans"
    if cis_col in df.columns and trans_col in df.columns:
        total_ct = df[cis_col] + df[trans_col]
        df["CisRatio"] = (df[cis_col] / total_ct).where(total_ct > 0, 0)
        df["TransRatio"] = (df[trans_col] / total_ct).where(total_ct > 0, 0)

    cis10kb_col = f"UniqMapped{mq}NoDup_Cis10kb"
    if cis10kb_col in df.columns and nodup_col in df.columns:
        df["LongRangeCisRate"] = (
            df[cis10kb_col] / df[nodup_col]
        ).where(df[nodup_col] > 0, 0)

    # ── Derived ratios (pairs-based) ──
    if "TotalPairs" in df.columns:
        if "Pairs_Cis" in df.columns and "Pairs_Trans" in df.columns:
            total_pt = df["Pairs_Cis"] + df["Pairs_Trans"]
            df["Pairs_CisRatio"] = (df["Pairs_Cis"] / total_pt).where(total_pt > 0, 0)
            df["Pairs_TransRatio"] = (df["Pairs_Trans"] / total_pt).where(total_pt > 0, 0)
        if "Pairs_Cis10kb" in df.columns:
            df["Pairs_LongRangeCisRate"] = (
                df["Pairs_Cis10kb"] / df["TotalPairs"]
            ).where(df["TotalPairs"] > 0, 0)

    # ── Reorder columns ──
    bam_cols = [
        "TotalFragments", "UniqMapped", f"UniqMapped{mq}",
        "Duplicates", f"UniqMapped{mq}NoDup",
        f"UniqMapped{mq}NoDup_Cis", f"UniqMapped{mq}NoDup_Trans",
        f"UniqMapped{mq}NoDup_Cis1kb", f"UniqMapped{mq}NoDup_Cis10kb",
        f"UniqMapped{mq}NoDup_Cis20kb",
        "MappingRate", f"{mq}_Rate", "DuplicationRate",
        "CisRatio", "TransRatio", "LongRangeCisRate",
    ]
    pairs_cols = [
        "TotalPairs", "Pairs_Cis", "Pairs_Trans",
        "Pairs_Cis1kb", "Pairs_Cis10kb", "Pairs_Cis20kb",
        "Pairs_CisRatio", "Pairs_TransRatio", "Pairs_LongRangeCisRate",
    ]
    # Pair type columns
    pt_cols = sorted([c for c in df.columns if c.startswith("PairType_")])

    ordered = [c for c in bam_cols + pairs_cols + pt_cols if c in df.columns]
    extra = [c for c in df.columns if c not in ordered]
    df = df[ordered + extra]

    return df


def sample_summary(df, mapq_cutoff=30):
    """Produce a one-row sample-level summary from per-cell DataFrame."""
    mq = f"MapQ{mapq_cutoff}"
    summary = {}

    summary["TotalBarcodes"] = len(df)

    # Sum up key counts
    for col in ["TotalFragments", "UniqMapped", f"UniqMapped{mq}",
                 f"UniqMapped{mq}NoDup",
                 f"UniqMapped{mq}NoDup_Cis", f"UniqMapped{mq}NoDup_Trans",
                 f"UniqMapped{mq}NoDup_Cis10kb",
                 "TotalPairs", "Pairs_Cis", "Pairs_Trans", "Pairs_Cis10kb"]:
        if col in df.columns:
            summary[f"Sum_{col}"] = int(df[col].sum())

    # Median per-cell stats
    for col in ["TotalFragments", "UniqMapped", f"UniqMapped{mq}NoDup",
                 "MappingRate", "DuplicationRate",
                 "CisRatio", "LongRangeCisRate",
                 "TotalPairs", "Pairs_CisRatio", "Pairs_LongRangeCisRate"]:
        if col in df.columns:
            summary[f"Median_{col}"] = float(df[col].median())
            summary[f"Mean_{col}"] = float(df[col].mean())

    return pd.Series(summary)


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Per-cell QC for Rupture Droplet Hi-C pipeline"
    )
    parser.add_argument("--bam", default=None,
                        help="CB-tagged BAM file (03.mapping/{sample}_{genome}.bam)")
    parser.add_argument("--pairs", default=None,
                        help="Deduplicated pairs file (03.mapping/{sample}_{genome}.sc.pairs.gz)")
    parser.add_argument("--output", required=True,
                        help="Output prefix. Produces {prefix}.per_cell.tsv and {prefix}.sample_summary.tsv")
    parser.add_argument("--mapq", type=int, default=30,
                        help="MAPQ cutoff (default: 30)")
    args = parser.parse_args()

    if args.bam is None and args.pairs is None:
        parser.error("At least one of --bam or --pairs is required")

    bam_counts = None
    pairs_counts = None

    if args.bam:
        bam_path = pathlib.Path(args.bam)
        if not bam_path.exists():
            print(f"WARNING: BAM not found: {bam_path}")
        else:
            bam_counts = qc_from_bam(bam_path, mapq_cutoff=args.mapq)

    if args.pairs:
        pairs_path = pathlib.Path(args.pairs)
        if not pairs_path.exists():
            print(f"WARNING: Pairs file not found: {pairs_path}")
        else:
            pairs_counts = qc_from_pairs(pairs_path)

    # Merge
    merged = merge_stats(bam_counts, pairs_counts)

    if not merged:
        print("ERROR: No data collected from any input.")
        return

    print(f"\nBuilding per-cell QC table for {len(merged):,} barcodes...")
    df = build_dataframe(merged, mapq_cutoff=args.mapq)


    # Save per-cell
    per_cell_out = f"{args.output}.per_cell.tsv"
    df.to_csv(per_cell_out, sep="\t")
    print(f"Saved: {per_cell_out}  ({df.shape[0]} barcodes x {df.shape[1]} columns)")

    # Save sample summary
    summary = sample_summary(df, mapq_cutoff=args.mapq)
    summary_out = f"{args.output}.sample_summary.tsv"
    summary.to_frame(name="Value").to_csv(summary_out, sep="\t")
    print(f"Saved: {summary_out}")

    # Print overview
    print("\n--- Per-cell overview ---")
    overview_cols = [
        "TotalFragments", f"UniqMappedMapQ{args.mapq}NoDup",
        "MappingRate", "DuplicationRate", "CisRatio", "LongRangeCisRate",
        "TotalPairs", "Pairs_CisRatio", "Pairs_LongRangeCisRate",
    ]
    for col in overview_cols:
        if col in df.columns:
            s = df[col].dropna()
            if len(s) > 0:
                print(f"  {col:45s}  median={s.median():>12.4f}  "
                      f"mean={s.mean():>12.4f}  "
                      f"min={s.min():>12.4f}  max={s.max():>12.4f}")


if __name__ == "__main__":
    main()