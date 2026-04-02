#!/usr/bin/env python
"""
Per-cell QC for snmC-seq2 pipeline (Bismark SE, PBAT R1 + non-directional R2).

Trimming was done with cutadapt (two-pass: adapter trim + 5' clip).
Cutadapt does not write a Trim Galore-style report, so trim stats are
collected from the Bismark report's input read count.

Parses:
  1. Bismark SE alignment reports  (05.align/{prefix}_{1,2}.clean_bismark_bt2_SE_report.txt)
  2. BAM summary txt               (05.align/{prefix}_{1,2}.summary.txt, from se_bam_summary.py)
  3. BisQC summary.txt             (05.align/{prefix}_{1,2}.rmdup.RG.summary.txt, if available)

Output:
  {output_dir}/{cell_id}.qc_stats.csv

Usage:
    python snmcseq2_qc_per_cell.py \\
        --cell_id SRR1234567 \\
        --project_dir /path/to/snmCseq2 \\
        --output_dir qc_stats \\
        --genome hg38        # or mm10
"""

import argparse
import pathlib
import re
from collections import OrderedDict

import pandas as pd


# ═══════════════════════════════════════════════════════════════════
# 1. Bismark SE alignment reports
# ═══════════════════════════════════════════════════════════════════

def parse_bismark_se_report(path):
    """Parse a Bismark SE alignment report."""
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Sequences analysed in total:" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["TotalReads"] = int(m.group(1))
            elif "unique best hit" in line and "Number of" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["UniqMapped"] = int(m.group(1))
            elif "Mapping efficiency:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["MappingRate"] = float(m.group(1))
            elif "Sequences with no alignments" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Unmapped"] = int(m.group(1))
            elif "not map uniquely" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Ambiguous"] = int(m.group(1))
            elif "Total number of C's analysed:" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["TotalC"] = int(m.group(1))
            elif "C methylated in CpG context:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["mCG_Rate"] = float(m.group(1))
            elif "C methylated in CHG context:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["mCHG_Rate"] = float(m.group(1))
            elif "C methylated in CHH context:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["mCHH_Rate"] = float(m.group(1))
    return stats


def collect_bismark_stats(project_dir, prefix):
    """
    Collect Bismark SE report for R1 (PBAT) and R2 (non-directional).

    Bismark on input {prefix}_1.clean.fq.gz (no --prefix flag) produces:
        05.align/{prefix}_1.clean_bismark_bt2_SE_report.txt
    """
    align_dir = project_dir / "05.align"
    stats = {}
    for read_num, label in [("1", "R1"), ("2", "R2")]:
        rpt = align_dir / f"{prefix}_{read_num}.clean_bismark_bt2_SE_report.txt"
        if not rpt.exists():
            print(f"  [WARN] Missing Bismark report: {rpt}")
        parsed = parse_bismark_se_report(rpt)
        for k, v in parsed.items():
            stats[f"Bismark_{label}_{k}"] = v
    return stats


# ═══════════════════════════════════════════════════════════════════
# 2. BAM summary txt (se_bam_summary.py output)
# ═══════════════════════════════════════════════════════════════════

def parse_bam_summary_txt(path):
    """
    Parse tab-separated summary txt from se_bam_summary.py or bam_summary_universal.py.

    Format (tab-separated):
        Key:    <int>    <float>%
    """
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            line = line.strip()
            if (not line
                    or line.startswith("#")
                    or line.upper().startswith("WARNING")
                    or line.upper().startswith("INFO")
                    or line.upper().startswith("NOTE")):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key = parts[0].rstrip(":")
            raw_val = parts[1].strip()
            try:
                stats[key] = int(raw_val)
            except ValueError:
                stats[key] = raw_val
            if len(parts) >= 3:
                pct_str = parts[2].strip().rstrip("%")
                try:
                    stats[key + "_Pct"] = float(pct_str)
                except ValueError:
                    pass
    return stats


def collect_bam_summary_stats(project_dir, prefix):
    """
    Collect se_bam_summary stats for R1 and R2.
    Note: run on raw (non-rmdup) BAMs, so DuplicateMapped reflects
    BAM-flagged duplicates in the raw alignment.
    """
    align_dir = project_dir / "05.align"
    stats = {}
    for read_num, label in [("1", "R1"), ("2", "R2")]:
        txt = align_dir / f"{prefix}_{read_num}.summary.txt"
        if not txt.exists():
            print(f"  [WARN] Missing BAM summary: {txt}")
        parsed = parse_bam_summary_txt(txt)
        for k, v in parsed.items():
            stats[f"BAM_{label}_{k}"] = v
    return stats


# ═══════════════════════════════════════════════════════════════════
# 3. Site counts (6plus2 bed) and trinuc conversion rates
# ═══════════════════════════════════════════════════════════════════

def _count_bed_rows(path):
    """Return number of data rows in a bed file."""
    if not path.exists():
        return None
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            if line.strip():
                n += 1
    return n


def _get_trinuc_rate(path, context):
    """Return methylation rate (%) for one trinuc context from a trinuc_methy.txt file."""
    if not path.exists():
        return None
    try:
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[0].rstrip(":") == context:
                    pct = parts[2].rstrip("%")
                    return float(pct) if pct not in ("NaN", "nan", "") else None
    except Exception as e:
        print(f"  [trinuc err {path.name}: {e}]")
    return None


def collect_site_and_trinuc_stats(project_dir, prefix):
    """
    Collect:
      - HCG_site_count : rows in R1+R2 CG.6plus2.bed  (covered CpG sites)
      - chrM_noncpg    : ACT rate on chrM (bisulfite conversion proxy)
      - chrM_endo      : ACG rate on chrM
      - chr21_noncpg   : ACT rate on chr21
      - chr21_endo     : ACG rate on chr21
    """
    align_dir = project_dir / "05.align"
    stats = {}

    # CG.6plus2.bed row count — sum R1 and R2
    total = 0
    for read_num in ("1", "2"):
        bed = align_dir / f"{prefix}_{read_num}.rmdup.RG.cpg.filtered.sort.CG.6plus2.bed"
        n = _count_bed_rows(bed)
        if n is not None:
            total += n
    stats["HCG_site_count"] = total if total > 0 else None

    # trinuc rates — average R1 and R2
    for chrom in ("chrM", "chr21"):
        for ctx, label in [("ACT", "noncpg"), ("ACG", "endo")]:
            vals = []
            for read_num in ("1", "2"):
                path = align_dir / f"{prefix}_{read_num}.rmdup.RG.trinuc_methy.{chrom}.txt"
                v = _get_trinuc_rate(path, ctx)
                if v is not None:
                    vals.append(v)
            stats[f"{chrom}_{label}"] = sum(vals) / len(vals) if vals else None

    return stats


# ═══════════════════════════════════════════════════════════════════
# 4. BisQC summary (optional, only if Bis-QC.pl was run)
# ═══════════════════════════════════════════════════════════════════

def parse_bisqc_summary(path):
    stats = {}
    if not path.exists():
        return stats
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                match = re.match(r'^(.+?)[\t:]\s*(.+)$', line)
                if match:
                    key = match.group(1).strip()
                    val = match.group(2).strip()
                    col = "BisQC_" + re.sub(r'[^A-Za-z0-9]+', '_', key).strip('_')
                    num_m = re.match(r'^([\d,.]+)', val.replace(',', ''))
                    if num_m:
                        try:
                            stats[col] = int(float(num_m.group(1)))
                        except ValueError:
                            stats[col] = val
                    else:
                        stats[col] = val
                    pct_m = re.search(r'([\d.]+)%', val)
                    if pct_m:
                        stats[col + "_Pct"] = float(pct_m.group(1))
    except Exception as e:
        print(f"  [bisqc err: {e}]")
    return stats


def collect_bisqc_stats(project_dir, prefix):
    align_dir = project_dir / "05.align"
    stats = {}
    for read_num, label in [("1", "R1"), ("2", "R2")]:
        path = align_dir / f"{prefix}_{read_num}.rmdup.RG.summary.txt"
        parsed = parse_bisqc_summary(path)
        for k, v in parsed.items():
            stats[f"{label}_{k}"] = v
    return stats


# ═══════════════════════════════════════════════════════════════════
# Derived combined metrics
# ═══════════════════════════════════════════════════════════════════

def add_derived_stats(stats):
    """Combine R1+R2 into summary metrics."""
    r1_total = stats.get("Bismark_R1_TotalReads") or 0
    r2_total = stats.get("Bismark_R2_TotalReads") or 0
    r1_uniq = stats.get("Bismark_R1_UniqMapped") or 0
    r2_uniq = stats.get("Bismark_R2_UniqMapped") or 0
    if r1_total + r2_total > 0:
        stats["Combined_TotalReads"] = r1_total + r2_total
        stats["Combined_UniqMapped"] = r1_uniq + r2_uniq
        stats["Combined_MappingRate"] = (r1_uniq + r2_uniq) / (r1_total + r2_total) * 100

    for metric, key in [("mCG_Rate", "Combined_mCG_Rate"),
                        ("mCHG_Rate", "Combined_mCHG_Rate"),
                        ("mCHH_Rate", "Combined_mCHH_Rate")]:
        vals = [stats.get(f"Bismark_R{n}_{metric}") for n in ["1", "2"]]
        vals = [v for v in vals if v is not None]
        if vals:
            stats[key] = sum(vals) / len(vals)

    # Total mapped reads from BAM summary (R1 + R2)
    r1_mapped = stats.get("BAM_R1_Mapped") or 0
    r2_mapped = stats.get("BAM_R2_Mapped") or 0
    if r1_mapped or r2_mapped:
        stats["BAM_Combined_Mapped"] = r1_mapped + r2_mapped

    r1_mapq = stats.get("BAM_R1_MappedMapQ30") or 0
    r2_mapq = stats.get("BAM_R2_MappedMapQ30") or 0
    if r1_mapq or r2_mapq:
        stats["BAM_Combined_MappedMapQ30"] = r1_mapq + r2_mapq

    return stats


# ═══════════════════════════════════════════════════════════════════
# Collect all for one cell
# ═══════════════════════════════════════════════════════════════════

def collect_cell_stats(prefix, project_dir, genome=None, mapq_cutoff=30):
    project_dir = pathlib.Path(project_dir)
    stats = OrderedDict({"CellID": prefix})
    if genome:
        stats["Genome"] = genome
    stats.update(collect_bismark_stats(project_dir, prefix))
    stats.update(collect_bam_summary_stats(project_dir, prefix))
    stats.update(collect_site_and_trinuc_stats(project_dir, prefix))
    stats.update(collect_bisqc_stats(project_dir, prefix))
    add_derived_stats(stats)
    return stats


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Per-cell QC for snmC-seq2 pipeline (SLURM array-friendly)"
    )
    parser.add_argument("--cell_id", required=True, help="Cell/sample ID")
    parser.add_argument("--project_dir", required=True, help="Root of snmCseq2 project dir")
    parser.add_argument("--output_dir", required=True, help="Dir for per-cell .qc_stats.csv")
    parser.add_argument("--genome", default=None,
                        help="Genome label (e.g. hg38 or mm10) stored as metadata column")
    parser.add_argument("--mapq", type=int, default=30, help="MAPQ cutoff (default: 30)")
    args = parser.parse_args()

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print(f"Processing cell: {args.cell_id}  genome={args.genome or 'unset'}")
    stats = collect_cell_stats(
        args.cell_id, args.project_dir,
        genome=args.genome,
        mapq_cutoff=args.mapq,
    )
    n = sum(1 for k, v in stats.items() if k not in ("CellID", "Genome") and v not in (None, "", 0))
    print(f"  Collected {n} non-zero metrics")

    out_path = pathlib.Path(args.output_dir) / f"{args.cell_id}.qc_stats.csv"
    pd.Series(stats).to_frame().T.set_index("CellID").to_csv(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
