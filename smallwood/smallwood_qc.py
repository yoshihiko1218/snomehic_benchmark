#!/usr/bin/env python
"""
QC Summary Collector for Smallwood scBS-seq Pipeline
====================================================

Collects per-cell QC from:
  1. Trim Galore reports         (03.trimmed_fastq)
  2. Bismark SE report - mm10    (05.align_mm10)
  3. BAM stats from rmdup BAM    (05.align_mm10, via pysam)
  4. Trinuc methylation - chrM   (05.align_mm10, conversion rate)
  5. Trinuc methylation - chr19  (05.align_mm10, non-CG)
  6. CpG counts from .6plus2.bed (05.align_mm10)
  7. Bis-QC summary.txt          (05.align_mm10)

Usage:
    python smallwood_qc.py \
        --project_dir /path/to/smallwood \
        --acc_list acc_list.txt \
        --output smallwood_qc_summary.csv \
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

    # flagstat
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

    # MAPQ iteration
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
#    Just extract ACT/ACG/GCT percentages directly from the file.
# ═══════════════════════════════════════════════════════════════════

def parse_trinuc_methy(path):
    """
    Parse trinuc_methy file and extract ACT (noncpg), ACG (endo), GCT (exo).
    Same logic as your extract_values script.
    """
    if not path.exists():
        return None, None, None

    try:
        methy_df = pd.read_csv(
            path, sep='\t', header=None,
            names=['trinuc', 'count', 'percent']
        )

        def get_percent(trinuc):
            vals = methy_df.loc[methy_df['trinuc'] == f"{trinuc}:", 'percent']
            return float(vals.values[0].strip('%')) if len(vals) > 0 else None

        noncpg = get_percent('ACT')
        endo = get_percent('ACG')
        exo = get_percent('GCT')
        return noncpg, endo, exo

    except Exception as e:
        print(f" [trinuc err: {e}]", end="")
        return None, None, None


def collect_conversion_stats(project_dir, prefix):
    align_dir = project_dir / "05.align_mm10"
    stats = {}

    # chrM
    noncpg, endo, exo = parse_trinuc_methy(
        align_dir / f"{prefix}.rmdup.RG.trinuc_methy.chrM.txt")
    stats["chrM_noncpg"] = noncpg
    stats["chrM_endo"] = endo
    stats["chrM_exo"] = exo

    # chr19
    noncpg, endo, exo = parse_trinuc_methy(
        align_dir / f"{prefix}.rmdup.RG.trinuc_methy.chr19.txt")
    stats["chr19_noncpg"] = noncpg
    stats["chr19_endo"] = endo
    stats["chr19_exo"] = exo

    return stats


# ═══════════════════════════════════════════════════════════════════
# 5. CpG counts from .6plus2.bed
# ═══════════════════════════════════════════════════════════════════

def collect_cpg_bed_stats(project_dir, prefix):
    """
    Parse .rmdup.RG.cpg.filtered.sort.CG.6plus2.bed

    Columns: chr start end . score strand methyl_pct read_count
    (skip track header line)
    """
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
# 7. Collect all for one cell
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
        description="Collect QC summary for Smallwood scBS-seq pipeline"
    )
    parser.add_argument("--project_dir", required=True)
    parser.add_argument("--acc_list", required=True)
    parser.add_argument("--output", default="smallwood_qc_summary.csv")
    parser.add_argument("--mapq", type=int, default=30)
    args = parser.parse_args()

    project_dir = pathlib.Path(args.project_dir).resolve()
    with open(args.acc_list) as f:
        cell_ids = [line.strip() for line in f if line.strip()]

    print(f"Project dir: {project_dir}")
    print(f"Cells: {len(cell_ids)}")

    records = []
    for prefix in cell_ids:
        print(f"  {prefix}...", end="")
        stats = collect_cell_stats(prefix, project_dir, mapq_cutoff=args.mapq)
        n = sum(1 for k, v in stats.items() if k != "CellID" and v not in (None, "", 0))
        print(f" {n} metrics")
        records.append(stats)

    df = pd.DataFrame(records).set_index("CellID")

    # Reorder
    priority_cols = [
        "Trim_R1_InputReads", "Trim_R1_TrimmedReads", "Trim_R1_TrimmedReads_Pct",
        "Trim_R2_InputReads", "Trim_R2_TrimmedReads", "Trim_R2_TrimmedReads_Pct",
        "Bismark_TotalReads", "Bismark_UniqMapped", "Bismark_MappingRate",
        "Bismark_Unmapped", "Bismark_Ambiguous",
        "Bismark_mCG_Rate", "Bismark_mCHG_Rate", "Bismark_mCHH_Rate",
        "Bismark_TotalC",
        "BAM_TotalAlignments", "BAM_PrimaryMapped", "BAM_Duplicates",
        f"BAM_MapQ{args.mapq}", f"BAM_MapQ{args.mapq}_Rate",
        "chrM_noncpg", "chrM_endo", "chrM_exo",
        "chr19_noncpg", "chr19_endo", "chr19_exo",
        "CpG_TotalSites", "CpG_TotalReads", "CpG_MeanMethylation",
        "CpG_MeanCoverage", "CpG_MethylatedSites", "CpG_MethylatedFrac",
    ]
    existing = [c for c in priority_cols if c in df.columns]
    extra = [c for c in df.columns if c not in priority_cols]
    df = df[existing + extra]

    df.to_csv(args.output)
    print(f"\nSaved: {args.output}")
    print(f"Shape: {df.shape[0]} cells x {df.shape[1]} columns")

    print("\n--- Key Metrics Overview ---")
    overview_cols = [
        "Trim_R1_InputReads", "Trim_R1_TrimmedReads_Pct",
        "Bismark_MappingRate", "Bismark_mCG_Rate",
        f"BAM_MapQ{args.mapq}", f"BAM_MapQ{args.mapq}_Rate",
        "chrM_noncpg", "chrM_endo", "chrM_exo",
        "chr19_noncpg", "chr19_endo", "chr19_exo",
        "CpG_TotalSites", "CpG_TotalReads", "CpG_MeanMethylation",
    ]
    for col in overview_cols:
        if col in df.columns:
            s = pd.to_numeric(df[col], errors="coerce").dropna()
            if len(s) > 0:
                print(f"  {col:40s}  median={s.median():>12.2f}  "
                      f"mean={s.mean():>12.2f}  "
                      f"range=[{s.min():.2f}, {s.max():.2f}]")


if __name__ == "__main__":
    main()
