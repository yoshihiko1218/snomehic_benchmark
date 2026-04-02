#!/usr/bin/env python
"""
Per-cell QC for Smallwood scBS-seq pipeline.
Designed to run as a SLURM array job (one task per cell).

Outputs: {output_dir}/{cell_id}.qc_stats.csv

Usage:
    python smallwood_qc_per_cell.py \
        --cell_id SRR1248444 \
        --project_dir /path/to/smallwood \
        --output_dir qc_stats \
        --mapq 30
"""

import argparse
import pathlib
import re
from collections import OrderedDict

import pandas as pd
import pysam


# ═══════════════════════════════════════════════════════════════════
# 1. Trim Galore reports
# ═══════════════════════════════════════════════════════════════════

def parse_trim_galore_report(path):
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
        rpt = trim_dir / f"{prefix}_{read_num}.fastq_trimming_report.txt"
        parsed = parse_trim_galore_report(rpt)
        for k, v in parsed.items():
            stats[f"Trim_{label}_{k}"] = v
    return stats


# ═══════════════════════════════════════════════════════════════════
# 2. Bismark SE report (mm10 only)
# ═══════════════════════════════════════════════════════════════════

def parse_bismark_se_report(path):
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Sequences analysed in total:" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_TotalReads"] = int(m.group(1))
            elif "unique best hit" in line and "Number of" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_UniqMapped"] = int(m.group(1))
            elif "Mapping efficiency:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["Bismark_MappingRate"] = float(m.group(1))
            elif "Sequences with no alignments" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_Unmapped"] = int(m.group(1))
            elif "not map uniquely" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_Ambiguous"] = int(m.group(1))
            elif line.startswith("CT/CT"):
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_OT"] = int(m.group(1))
            elif line.startswith("CT/GA"):
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_OB"] = int(m.group(1))
            elif line.startswith("GA/CT"):
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_CTOT"] = int(m.group(1))
            elif line.startswith("GA/GA"):
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_CTOB"] = int(m.group(1))
            elif "Total number of C's analysed:" in line:
                m = re.search(r'(\d+)', line.split(":")[-1])
                if m:
                    stats["Bismark_TotalC"] = int(m.group(1))
            elif "C methylated in CpG context:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["Bismark_mCG_Rate"] = float(m.group(1))
            elif "C methylated in CHG context:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["Bismark_mCHG_Rate"] = float(m.group(1))
            elif "C methylated in CHH context:" in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    stats["Bismark_mCHH_Rate"] = float(m.group(1))
    return stats


def collect_bismark_mm10_se(project_dir, prefix):
    rpt = project_dir / "05.align_mm10" / f"{prefix}.SE.input_bismark_bt2_SE_report.txt"
    return parse_bismark_se_report(rpt)


# ═══════════════════════════════════════════════════════════════════
# 3. BAM stats (from rmdup.RG.bam via pysam)
# ═══════════════════════════════════════════════════════════════════

def collect_bam_stats(project_dir, prefix, mapq_cutoff=30):
    bam_path = project_dir / "05.align_mm10" / f"{prefix}.rmdup.RG.bam"
    stats = {}
    if not bam_path.exists():
        return stats

    try:
        flagstat_text = pysam.flagstat(str(bam_path))
        for line in flagstat_text.strip().split("\n"):
            count = int(line.split(" + ")[0])
            if "in total" in line:
                stats["BAM_TotalAlignments"] = count
            elif "secondary" in line:
                stats["BAM_Secondary"] = count
            elif "supplementary" in line:
                stats["BAM_Supplementary"] = count
            elif "duplicates" in line and "primary duplicates" not in line:
                stats["BAM_Duplicates"] = count
            elif "primary mapped" in line:
                stats["BAM_PrimaryMapped"] = count
            elif "mapped" in line and "mate mapped" not in line and "primary mapped" not in line:
                stats["BAM_Mapped"] = count
    except Exception as e:
        print(f" [flagstat err: {e}]", end="")

    try:
        mapped_primary = 0
        mapq_pass = 0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                mapped_primary += 1
                if read.mapping_quality >= mapq_cutoff:
                    mapq_pass += 1
        stats["BAM_PrimaryMapped_iter"] = mapped_primary
        stats[f"BAM_MapQ{mapq_cutoff}"] = mapq_pass
        if mapped_primary > 0:
            stats[f"BAM_MapQ{mapq_cutoff}_Rate"] = mapq_pass / mapped_primary * 100
    except Exception as e:
        print(f" [bam iter err: {e}]", end="")

    return stats


# ═══════════════════════════════════════════════════════════════════
# 4. Trinuc methylation (chrM + chr19)
# ═══════════════════════════════════════════════════════════════════

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
        print(f" [trinuc err: {e}]", end="")
    return None


def collect_conversion_stats(project_dir, prefix):
    """
    Extract noncpg / endo / exo rates from chrM and chr19 trinuc files.

      noncpg  — ACT rate (non-CG, bisulfite conversion proxy)
      endo    — ACG rate (endogenous CpG methylation proxy)
      exo     — GCT rate (GCH, NOMe accessibility proxy)
    """
    align_dir = project_dir / "05.align_mm10"
    stats = {}
    for chrom in ("chrM", "chr19"):
        path = align_dir / f"{prefix}.rmdup.RG.trinuc_methy.{chrom}.txt"
        for ctx, label in [("ACT", "noncpg"), ("ACG", "endo"), ("GCT", "exo")]:
            stats[f"{chrom}_{label}"] = _get_trinuc_rate(path, ctx)
    return stats


# ═══════════════════════════════════════════════════════════════════
# 5. CpG counts from .6plus2.bed
# ═══════════════════════════════════════════════════════════════════

def collect_cpg_bed_stats(project_dir, prefix):
    bed_path = (project_dir / "05.align_mm10" /
                f"{prefix}.rmdup.RG.cpg.filtered.sort.CG.6plus2.bed")
    stats = {}
    if not bed_path.exists():
        return stats

    try:
        total_sites = 0
        total_reads = 0
        total_meth_pct = 0.0
        methylated_sites = 0

        with open(bed_path) as f:
            for line in f:
                if line.startswith("track ") or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                meth_pct = float(parts[6])
                read_count = int(parts[7])

                total_sites += 1
                total_reads += read_count
                total_meth_pct += meth_pct
                if meth_pct > 0:
                    methylated_sites += 1

        stats["CpG_TotalSites"] = total_sites
        stats["CpG_TotalReads"] = total_reads
        stats["CpG_MethylatedSites"] = methylated_sites
        if total_sites > 0:
            stats["CpG_MeanMethylation"] = total_meth_pct / total_sites
            stats["CpG_MeanCoverage"] = total_reads / total_sites
            stats["CpG_MethylatedFrac"] = methylated_sites / total_sites * 100
    except Exception as e:
        print(f" [bed err: {e}]", end="")

    return stats


# ═══════════════════════════════════════════════════════════════════
# 6. Bis-QC summary.txt
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
        print(f" [summary err: {e}]", end="")
    return stats


# ═══════════════════════════════════════════════════════════════════
# Collect all for one cell
# ═══════════════════════════════════════════════════════════════════

def collect_cell_stats(prefix, project_dir, mapq_cutoff=30):
    project_dir = pathlib.Path(project_dir)
    stats = OrderedDict({"CellID": prefix})
    stats.update(collect_trim_stats(project_dir, prefix))
    stats.update(collect_bismark_mm10_se(project_dir, prefix))
    stats.update(collect_bam_stats(project_dir, prefix, mapq_cutoff))
    stats.update(collect_conversion_stats(project_dir, prefix))
    stats.update(collect_cpg_bed_stats(project_dir, prefix))
    stats.update(parse_bisqc_summary(
        project_dir / "05.align_mm10" / f"{prefix}.summary.txt"))
    return stats


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Per-cell QC for Smallwood scBS-seq pipeline (SLURM-friendly)"
    )
    parser.add_argument("--cell_id", required=True)
    parser.add_argument("--project_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--mapq", type=int, default=30)
    args = parser.parse_args()

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print(f"Processing cell: {args.cell_id}")
    stats = collect_cell_stats(args.cell_id, args.project_dir, mapq_cutoff=args.mapq)

    out_path = pathlib.Path(args.output_dir) / f"{args.cell_id}.qc_stats.csv"
    pd.Series(stats).to_frame().T.set_index("CellID").to_csv(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
