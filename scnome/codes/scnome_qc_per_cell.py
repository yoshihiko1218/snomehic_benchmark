#!/usr/bin/env python
"""
Per-cell QC for scNOMe-seq pipeline (hg38, Trim Galore + Bismark SE).

R1 and R2 are trimmed and aligned independently as single-end reads.

Parses:
  1. Trim Galore trimming reports  (03.trimmed_fastq/{prefix}_{1,2}.fastq_trimming_report.txt)
  2. Bismark SE alignment reports  (04.alignment/{prefix}.{prefix}_{1,2}_trimmed_bismark_bt2_SE_report.txt)
  3. BAM summary txt               (04.alignment/{prefix}_{1,2}.summary.txt, from bam_summary_universal.py)
  4. BisQC summary.txt             (04.alignment/{prefix}_{1,2}.rmdup.RG.summary.txt, if available)

Output:
  {output_dir}/{cell_id}.qc_stats.csv

Usage:
    python scnome_qc_per_cell.py \\
        --cell_id SRR12345678 \\
        --project_dir /path/to/scnome \\
        --output_dir qc_stats \\
        --mapq 30
"""

import argparse
import pathlib
import re
from collections import OrderedDict

import pandas as pd


# ═══════════════════════════════════════════════════════════════════
# 1. Trim Galore reports
# ═══════════════════════════════════════════════════════════════════

def parse_trim_galore_report(path):
    """Parse a Trim Galore trimming report, return dict of stats."""
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Total reads processed:" in line:
                m = re.search(r'([\d,]+)', line.split(":")[-1])
                if m:
                    stats["InputReads"] = int(m.group(1).replace(',', ''))
            elif "Reads with adapters:" in line:
                m = re.search(r'([\d,]+)', line.split(":")[-1])
                if m:
                    stats["ReadsWithAdapters"] = int(m.group(1).replace(',', ''))
                pct = re.search(r'\(([\d.]+)%\)', line)
                if pct:
                    stats["ReadsWithAdapters_Pct"] = float(pct.group(1))
            elif "Quality-trimmed:" in line:
                m = re.search(r'([\d,]+)\s+bp', line)
                if m:
                    stats["QualityTrimmedBP"] = int(m.group(1).replace(',', ''))
            elif "Reads written (passing filters):" in line:
                m = re.search(r'([\d,]+)', line.split(":")[-1])
                if m:
                    stats["TrimmedReads"] = int(m.group(1).replace(',', ''))
                pct = re.search(r'\(([\d.]+)%\)', line)
                if pct:
                    stats["TrimmedReads_Pct"] = float(pct.group(1))
            elif "Total written (filtered):" in line:
                m = re.search(r'([\d,]+)\s+bp', line)
                if m:
                    stats["TrimmedBP"] = int(m.group(1).replace(',', ''))
    return stats


def collect_trim_stats(project_dir, prefix):
    trim_dir = project_dir / "03.trimmed_fastq"
    stats = {}
    for read_num, label in [("1", "R1"), ("2", "R2")]:
        # Trim Galore names the report after the input file. Cells re-trimmed
        # from .fq.gz inputs get "<p>_<n>.fq.gz_trimming_report.txt"; older runs
        # from .fastq inputs get "<p>_<n>.fastq_trimming_report.txt". Try both.
        candidates = [
            trim_dir / f"{prefix}_{read_num}.fq.gz_trimming_report.txt",
            trim_dir / f"{prefix}_{read_num}.fastq_trimming_report.txt",
        ]
        rpt = next((c for c in candidates if c.exists()), candidates[-1])
        if not rpt.exists():
            print(f"  [WARN] Missing trim report: {rpt}")
        parsed = parse_trim_galore_report(rpt)
        for k, v in parsed.items():
            stats[f"Trim_{label}_{k}"] = v
    return stats


# ═══════════════════════════════════════════════════════════════════
# 2. Bismark SE alignment reports
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
    Collect Bismark SE report for R1 and R2.

    Bismark with --prefix PREFIX on input {prefix}_{N}_trimmed.fq.gz produces:
        04.alignment/{prefix}.{prefix}_{N}_trimmed_bismark_bt2_SE_report.txt
    """
    align_dir = project_dir / "04.alignment"
    stats = {}
    for read_num, label in [("1", "R1"), ("2", "R2")]:
        rpt = align_dir / f"{prefix}.{prefix}_{read_num}_trimmed_bismark_bt2_SE_report.txt"
        if not rpt.exists():
            print(f"  [WARN] Missing Bismark report: {rpt}")
        parsed = parse_bismark_se_report(rpt)
        for k, v in parsed.items():
            stats[f"Bismark_{label}_{k}"] = v
    return stats


# ═══════════════════════════════════════════════════════════════════
# 3. BAM summary txt (bam_summary_universal.py output)
# ═══════════════════════════════════════════════════════════════════

def parse_bam_summary_txt(path):
    """
    Parse tab-separated summary txt from bam_summary_universal.py or se_bam_summary.py.

    Format:
        Key:    <int>    <float>%       (count + optional percentage)
        Key:    <str>                   (non-numeric values)
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
    align_dir = project_dir / "04.alignment"
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
# 4. Site counts (6plus2 bed) and trinuc conversion rates
# ═══════════════════════════════════════════════════════════════════

def _count_bed_rows(path):
    """Return number of data rows in a bed file (ignores header/comment lines)."""
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
    """
    Extract the methylation rate (%) for a single trinucleotide context from a
    trinuc_methy.txt file.  Returns None if the file or context is absent.

    File format:  ACG:    <count>    <rate>%
    """
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
      - HCG_site_count : rows in R1+R2 HCG.6plus2.bed  (covered CpG sites)
      - GCH_site_count : rows in R1+R2 GCH.6plus2.bed  (NOMe accessibility sites)
      - chrM_noncpg    : ACT rate on chrM (bisulfite conversion proxy)
      - chrM_endo      : ACG rate on chrM
      - chrM_exo       : GCT rate on chrM (GCH, NOMe)
      - chr21_noncpg   : ACT rate on chr21
      - chr21_endo     : ACG rate on chr21
      - chr21_exo      : GCT rate on chr21 (GCH, NOMe)
    """
    align_dir = project_dir / "04.alignment"
    stats = {}

    # 6plus2 bed row counts — sum R1 and R2
    for context_label, bed_suffix in [("HCG", "HCG"), ("GCH", "GCH")]:
        total = 0
        for read_num in ("1", "2"):
            bed = align_dir / (
                f"{prefix}_{read_num}.rmdup.RG.cytosine.filtered.sort.{bed_suffix}.6plus2.bed"
            )
            n = _count_bed_rows(bed)
            if n is not None:
                total += n
        stats[f"{context_label}_site_count"] = total if total > 0 else None

    # trinuc rates — average R1 and R2 for each chrom
    for chrom in ("chrM", "chr21"):
        for ctx, label in [("ACT", "noncpg"), ("ACG", "endo"), ("GCT", "exo")]:
            vals = []
            for read_num in ("1", "2"):
                path = align_dir / f"{prefix}_{read_num}.rmdup.RG.trinuc_methy.{chrom}.txt"
                v = _get_trinuc_rate(path, ctx)
                if v is not None:
                    vals.append(v)
            stats[f"{chrom}_{label}"] = sum(vals) / len(vals) if vals else None

    return stats


# ═══════════════════════════════════════════════════════════════════
# 5. BisQC summary (optional, only if Bis-QC.pl was run)
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
    align_dir = project_dir / "04.alignment"
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
    """Add combined R1+R2 summary metrics."""
    r1_input = stats.get("Trim_R1_InputReads") or 0
    r2_input = stats.get("Trim_R2_InputReads") or 0
    if r1_input or r2_input:
        stats["TotalInputReads"] = r1_input + r2_input

    r1_total = stats.get("Bismark_R1_TotalReads") or 0
    r2_total = stats.get("Bismark_R2_TotalReads") or 0
    r1_uniq = stats.get("Bismark_R1_UniqMapped") or 0
    r2_uniq = stats.get("Bismark_R2_UniqMapped") or 0
    if r1_total + r2_total > 0:
        stats["Combined_TotalReads"] = r1_total + r2_total
        stats["Combined_UniqMapped"] = r1_uniq + r2_uniq
        stats["Combined_MappingRate"] = (r1_uniq + r2_uniq) / (r1_total + r2_total) * 100

    # Average mCG across R1 and R2
    r1_mcg = stats.get("Bismark_R1_mCG_Rate")
    r2_mcg = stats.get("Bismark_R2_mCG_Rate")
    vals = [v for v in [r1_mcg, r2_mcg] if v is not None]
    if vals:
        stats["Combined_mCG_Rate"] = sum(vals) / len(vals)

    r1_mchg = stats.get("Bismark_R1_mCHG_Rate")
    r2_mchg = stats.get("Bismark_R2_mCHG_Rate")
    vals = [v for v in [r1_mchg, r2_mchg] if v is not None]
    if vals:
        stats["Combined_mCHG_Rate"] = sum(vals) / len(vals)

    r1_mchh = stats.get("Bismark_R1_mCHH_Rate")
    r2_mchh = stats.get("Bismark_R2_mCHH_Rate")
    vals = [v for v in [r1_mchh, r2_mchh] if v is not None]
    if vals:
        stats["Combined_mCHH_Rate"] = sum(vals) / len(vals)

    return stats


# ═══════════════════════════════════════════════════════════════════
# Collect all for one cell
# ═══════════════════════════════════════════════════════════════════

def collect_cell_stats(prefix, project_dir, mapq_cutoff=30):
    project_dir = pathlib.Path(project_dir)
    stats = OrderedDict({"CellID": prefix})
    stats.update(collect_trim_stats(project_dir, prefix))
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
        description="Per-cell QC for scNOMe pipeline (SLURM array-friendly)"
    )
    parser.add_argument("--cell_id", required=True, help="Cell/sample ID")
    parser.add_argument("--project_dir", required=True, help="Root of scnome project dir")
    parser.add_argument("--output_dir", required=True, help="Dir for per-cell .qc_stats.csv")
    parser.add_argument("--mapq", type=int, default=30, help="MAPQ cutoff (default: 30)")
    args = parser.parse_args()

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print(f"Processing cell: {args.cell_id}")
    stats = collect_cell_stats(args.cell_id, args.project_dir, mapq_cutoff=args.mapq)
    n = sum(1 for k, v in stats.items() if k != "CellID" and v not in (None, "", 0))
    print(f"  Collected {n} non-zero metrics")

    out_path = pathlib.Path(args.output_dir) / f"{args.cell_id}.qc_stats.csv"
    pd.Series(stats).to_frame().T.set_index("CellID").to_csv(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
