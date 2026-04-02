#!/usr/bin/env python
"""
Unified QC comparison across all 9 benchmark methods.

Reads per-method QC summary CSVs, normalises them to canonical column names,
and produces three PDF figure sets:
  1. Base alignment metrics  — all methods
  2. Methylation metrics     — bisulfite methods only (yap, smallwood, scnome, snmcseq2)
  3. HiC contact metrics     — HiC/3C methods only  (yap m3c/scnomehic, droplethic)

Config file (CSV or TSV, with header):
    path,label,type
    /path/to/methylhic/alignment/stats/MappingSummary.csv.gz,MethylHiC,yap
    /path/to/methylhic_new/alignment/stats/MappingSummary.csv.gz,MethylHiC-new,yap
    /path/to/snmCseq3/alignment/stats/MappingSummary.csv.gz,snmC-seq3,yap
    /path/to/scnomehic/alignment/stats/MappingSummary.csv.gz,scNOME-HiC,yap
    /path/to/smallwood/smallwood_qc_summary.csv,Smallwood,smallwood
    /path/to/scnome/scnome_qc_summary.csv,scNOMe,scnome
    /path/to/snmCseq2/snmcseq2_qc_summary.csv,snmC-seq2,snmcseq2
    /path/to/nagano/qc_summary.csv,Nagano,nagano
    /path/to/droplethic/qc_stats/SRR27586278_hg38.per_cell.tsv,DropletHiC,droplethic

Supported type values:
    yap        — YAP/m3C MappingSummary.csv.gz
    smallwood  — Smallwood scBS-seq (smallwood_qc_summary.csv)
    scnome     — scNOMe-seq (scnome_qc_summary.csv)
    snmcseq2   — snmC-seq2  (snmcseq2_qc_summary.csv)
    nagano     — Nagano scHiC (qc_summary.csv)
    droplethic — Droplet Hi-C (per_cell.tsv)

Usage:
    python compare_qc_all.py \\
        --config datasets.csv \\
        --output all_methods_qc

Outputs:
    all_methods_qc.base_alignment.pdf
    all_methods_qc.methylation.pdf
    all_methods_qc.hic_contacts.pdf
    all_methods_qc.harmonized.csv
    all_methods_qc.summary.csv
"""

import argparse
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════════════
# Canonical metric groups
# ═══════════════════════════════════════════════════════════════════

BASE_ALIGNMENT = [
    ("InputReads",      "Input Reads"),
    ("MapQ30Rate",      "MapQ≥30 Rate (%)"),
    ("DupRate",         "Duplication Rate (%)"),
    ("FinalMappedReads","Final Mapped Reads"),
]

METHYLATION = [
    ("mCG_Rate",        "mCG Rate (%)"),
    ("mCH_Rate",        "mCH Rate (%)"),
    ("mCCC_Rate",       "mCCC Rate (%)"),      # YAP only
    ("LambdaConvProxy", "BS Conv. Proxy (%)"), # lambda or chrM noncpg
    ("NonCG_rate",      "Non-CG Meth. (chrM+chr21 ACT %)"),
    ("GCH_rate",        "GCH Meth. (NOMe, chr21 GCT %)"),
]

HIC_CONTACTS = [
    ("TotalContacts",   "Total Contacts"),
    ("CisLongRatio_Pct","Cis Long Ratio (%)"),
    ("TransRatio_Pct",  "Trans Ratio (%)"),
    ("LongRangeCisRate_Pct", "Long-Range Cis >10 kb (%)"),
]

# Which method types carry each modality
METHYL_TYPES = {"yap", "smallwood", "scnome", "snmcseq2"}
HIC_TYPES    = {"yap", "droplethic"}

# Methods that have lambda/conversion proxy
LAMBDA_TYPES = {"yap", "smallwood", "scnome", "snmcseq2"}


# ═══════════════════════════════════════════════════════════════════
# Harmonization: normalise each method's native columns → canonical
# ═══════════════════════════════════════════════════════════════════

def _safe(df, col, fallback=None):
    """Return df[col] as numeric series, or a zero series if absent."""
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce")
    return pd.Series(fallback if fallback is not None else np.nan, index=df.index)


def harmonize_yap(df):
    """
    YAP MappingSummary columns (rates already 0-100, fracs 0-1):
        CellInputReadPairs, R1InputReads, R2InputReads
        R1MappingRateMapQ30, R2MappingRateMapQ30  [0-100]
        R1DuplicationRate, R2DuplicationRate       [0-100]
        LambdaCYFrac, mCGFrac, mCHFrac, mCCCFrac  [0-1]
        TotalContacts, CisLongRatio, TransRatio    [counts / 0-1]
        FinalmCReads
    """
    h = pd.DataFrame(index=df.index)
    h["InputReads"]       = (_safe(df, "R1InputReads") + _safe(df, "R2InputReads"))
    h["MapQ30Rate"]       = (_safe(df, "R1MappingRateMapQ30") + _safe(df, "R2MappingRateMapQ30")) / 2
    h["DupRate"]          = (_safe(df, "R1DuplicationRate") + _safe(df, "R2DuplicationRate")) / 2
    h["FinalMappedReads"] = _safe(df, "FinalmCReads")
    h["mCG_Rate"]         = _safe(df, "mCGFrac") * 100
    h["mCH_Rate"]         = _safe(df, "mCHFrac") * 100
    h["mCCC_Rate"]        = _safe(df, "mCCCFrac") * 100
    h["LambdaConvProxy"]  = _safe(df, "LambdaCYFrac") * 100
    h["TotalContacts"]    = _safe(df, "TotalContacts")
    h["CisLongRatio_Pct"] = _safe(df, "CisLongRatio") * 100
    h["TransRatio_Pct"]   = _safe(df, "TransRatio") * 100
    return h


def harmonize_smallwood(df):
    """
    Smallwood columns:
        Trim_R1_InputReads, Trim_R2_InputReads
        BAM_MapQ30_Rate [0-100], BAM_MapQ30 (count)
        BAM_Duplicates, BAM_PrimaryMapped
        Bismark_mCG_Rate, Bismark_mCHH_Rate, Bismark_mCHG_Rate [0-100]
        chrM_noncpg (non-CpG methylation on chrM = bisulfite conversion proxy)
    """
    h = pd.DataFrame(index=df.index)
    h["InputReads"]       = _safe(df, "Trim_R1_InputReads") + _safe(df, "Trim_R2_InputReads")
    h["MapQ30Rate"]       = _safe(df, "BAM_MapQ30_Rate")
    dup  = _safe(df, "BAM_Duplicates")
    prim = _safe(df, "BAM_PrimaryMapped").replace(0, np.nan)
    h["DupRate"]          = dup / prim * 100
    h["FinalMappedReads"] = _safe(df, "BAM_MapQ30")
    h["mCG_Rate"]         = _safe(df, "Bismark_mCG_Rate")
    # Average CHG + CHH as a single mCH proxy, or use CHH only (more common)
    chh = _safe(df, "Bismark_mCHH_Rate")
    chg = _safe(df, "Bismark_mCHG_Rate")
    h["mCH_Rate"]         = (chh + chg) / 2
    # chrM noncpg methylation = non-converted cytosines on chrM (lower = better conversion)
    h["LambdaConvProxy"]  = _safe(df, "chrM_noncpg")

    # CpG site count from 6plus2 bed (already computed in per-cell script)
    h["HCG_site_count"] = _safe(df, "CpG_TotalSites")

    return h


def harmonize_scnome(df):
    """
    scNOMe QC columns (new, two-read SE):
        TotalInputReads (= Trim R1 + R2 input)
        BAM_R1_UniqMappedMapQ30Reads, BAM_R2_UniqMappedMapQ30Reads
        BAM_R1_TotalReads (= total units in BAM)
        BAM_R1_DuplicatePrimaryMapped, BAM_R1_PrimaryMappedReads
        Bismark_R1_mCG_Rate, Bismark_R2_mCG_Rate [0-100 from Bismark]
        Bismark_R1_mCHH_Rate, Bismark_R2_mCHH_Rate
        BAM_R1_UniqMappedMapQ30ReadsToLambda_Pct  (% of MAPQ30 reads → lambda)
    """
    h = pd.DataFrame(index=df.index)
    h["InputReads"] = _safe(df, "TotalInputReads")
    if h["InputReads"].isna().all():
        h["InputReads"] = (_safe(df, "Bismark_R1_TotalReads")
                           + _safe(df, "Bismark_R2_TotalReads"))

    # MAPQ30 reads (R1 + R2) / total BAM reads
    mapq_r1 = _safe(df, "BAM_R1_UniqMappedMapQ30Reads")
    mapq_r2 = _safe(df, "BAM_R2_UniqMappedMapQ30Reads")
    bam_r1  = _safe(df, "BAM_R1_TotalReads").replace(0, np.nan)
    bam_r2  = _safe(df, "BAM_R2_TotalReads").replace(0, np.nan)
    h["MapQ30Rate"] = ((mapq_r1 / bam_r1 + mapq_r2 / bam_r2) / 2 * 100)

    dup_r1  = _safe(df, "BAM_R1_DuplicatePrimaryMapped")
    prim_r1 = _safe(df, "BAM_R1_PrimaryMappedReads").replace(0, np.nan)
    dup_r2  = _safe(df, "BAM_R2_DuplicatePrimaryMapped")
    prim_r2 = _safe(df, "BAM_R2_PrimaryMappedReads").replace(0, np.nan)
    h["DupRate"] = ((dup_r1 / prim_r1 + dup_r2 / prim_r2) / 2 * 100)

    h["FinalMappedReads"] = mapq_r1 + mapq_r2

    h["mCG_Rate"] = (_safe(df, "Bismark_R1_mCG_Rate") + _safe(df, "Bismark_R2_mCG_Rate")) / 2
    chh_r1 = _safe(df, "Bismark_R1_mCHH_Rate")
    chh_r2 = _safe(df, "Bismark_R2_mCHH_Rate")
    chg_r1 = _safe(df, "Bismark_R1_mCHG_Rate")
    chg_r2 = _safe(df, "Bismark_R2_mCHG_Rate")
    h["mCH_Rate"] = ((chh_r1 + chh_r2 + chg_r1 + chg_r2) / 4)

    # Lambda fraction from BAM summary (% of MAPQ30 reads mapping to lambda)
    h["LambdaConvProxy"] = _safe(df, "BAM_R1_UniqMappedMapQ30ReadsToLambda_Pct")

    # Site counts from 6plus2 bed (rows = covered sites, populated by per-cell script)
    h["HCG_site_count"] = _safe(df, "HCG_site_count")
    h["GCH_site_count"] = _safe(df, "GCH_site_count")

    return h


def harmonize_snmcseq2(df):
    """
    snmC-seq2 QC columns (new, two-read SE, cutadapt trim):
        Combined_TotalReads, Combined_MappingRate
        BAM_R1_MappedMapQ30_Pct [0-100], BAM_R2_MappedMapQ30_Pct
        BAM_Combined_MappedMapQ30
        BAM_R1_DuplicateMapped_Pct [0-100]
        Combined_mCG_Rate, Combined_mCHH_Rate, Combined_mCHG_Rate [0-100]
        BAM_R1_MappedMapQ30ToLambda_Pct
    """
    h = pd.DataFrame(index=df.index)
    h["InputReads"]       = _safe(df, "Combined_TotalReads")

    mapq_r1 = _safe(df, "BAM_R1_MappedMapQ30_Pct")
    mapq_r2 = _safe(df, "BAM_R2_MappedMapQ30_Pct")
    h["MapQ30Rate"]       = (mapq_r1 + mapq_r2) / 2

    dup_r1  = _safe(df, "BAM_R1_DuplicateMapped_Pct")
    dup_r2  = _safe(df, "BAM_R2_DuplicateMapped_Pct")
    h["DupRate"]          = (dup_r1 + dup_r2) / 2

    h["FinalMappedReads"] = _safe(df, "BAM_Combined_MappedMapQ30")
    h["mCG_Rate"]         = _safe(df, "Combined_mCG_Rate")
    chh = _safe(df, "Combined_mCHH_Rate")
    chg = _safe(df, "Combined_mCHG_Rate")
    h["mCH_Rate"]         = (chh + chg) / 2
    h["LambdaConvProxy"]  = (
        _safe(df, "BAM_R1_MappedMapQ30ToLambda_Pct") +
        _safe(df, "BAM_R2_MappedMapQ30ToLambda_Pct")
    ) / 2

    # Site counts from 6plus2 bed
    h["HCG_site_count"] = _safe(df, "HCG_site_count")

    return h


def harmonize_nagano(df):
    """
    Nagano scHiC (fastp + Bowtie2 PE):
        InputReadPairs, InputReads_R1, InputReads_R2
        MappingRateMapQ30 [0-100]
        DuplicationRate   [0-100]
        DeduppedReads, FinalReadPairs
    Note: no methylation, no contact ratios collected.
    """
    h = pd.DataFrame(index=df.index)
    inp = _safe(df, "InputReads_R1") + _safe(df, "InputReads_R2")
    if inp.isna().all():
        inp = _safe(df, "InputReadPairs") * 2
    h["InputReads"]       = inp
    h["MapQ30Rate"]       = _safe(df, "MappingRateMapQ30")
    h["DupRate"]          = _safe(df, "DuplicationRate")
    h["FinalMappedReads"] = _safe(df, "DeduppedReads")
    return h


def harmonize_droplethic(df):
    """
    Droplet HiC per-cell TSV (CB-indexed):
        TotalFragments
        MapQ30_Rate      [0-1]
        DuplicationRate  [0-1]
        UniqMappedMapQ30NoDup
        CisRatio, TransRatio, LongRangeCisRate  [0-1]
    """
    h = pd.DataFrame(index=df.index)
    h["InputReads"]           = _safe(df, "TotalFragments") * 2
    h["MapQ30Rate"]           = _safe(df, "MapQ30_Rate") * 100
    h["DupRate"]              = _safe(df, "DuplicationRate") * 100
    h["FinalMappedReads"]     = _safe(df, "UniqMappedMapQ30NoDup")
    # Use CisRatio as CisLong proxy (all valid pairs) and Cis10kb rate as long-range
    h["TotalContacts"]        = _safe(df, "UniqMappedMapQ30NoDup")
    h["CisLongRatio_Pct"]     = _safe(df, "CisRatio") * 100
    h["TransRatio_Pct"]       = _safe(df, "TransRatio") * 100
    h["LongRangeCisRate_Pct"] = _safe(df, "LongRangeCisRate") * 100
    return h


HARMONIZERS = {
    "yap":        harmonize_yap,
    "smallwood":  harmonize_smallwood,
    "scnome":     harmonize_scnome,
    "snmcseq2":   harmonize_snmcseq2,
    "nagano":     harmonize_nagano,
    "droplethic": harmonize_droplethic,
}


# ═══════════════════════════════════════════════════════════════════
# Load + harmonize each dataset
# ═══════════════════════════════════════════════════════════════════

def load_dataset(path, label, dtype):
    path = pathlib.Path(path)
    if not path.exists():
        print(f"  {label:20s}  [SKIP — file not found: {path}]")
        return None
    sep = "\t" if path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(path, index_col=0, sep=sep)
    df.insert(0, "Dataset", label)
    df.insert(1, "Type", dtype)
    print(f"  {label:20s}  {len(df):5d} cells  [{dtype}]  {path.name}")
    return df


def harmonize_dataset(df, dtype, label):
    if dtype not in HARMONIZERS:
        raise ValueError(f"Unknown type '{dtype}'. Valid: {list(HARMONIZERS)}")
    h = HARMONIZERS[dtype](df)
    h.insert(0, "Dataset", label)
    h.insert(1, "Type", dtype)
    h["CellID"] = df.index.astype(str)
    return h


# ═══════════════════════════════════════════════════════════════════
# Trinuc: load pre-computed noncpg / endo / exo from summary/trinuc/
# ═══════════════════════════════════════════════════════════════════

# Mapping from dataset type → trinuc file basename (relative to config dir's trinuc/ subfolder)
TRINUC_FILE_MAP = {
    "scnome":    "trinuc/scnome.chr21.txt",
    "snmcseq2":  "trinuc/snmCseq2.chr21.txt",
    "smallwood": "trinuc/smallwood.chr21.txt",
}

# Types whose sample IDs have a read-number suffix (_1 or _2) that must be stripped
TRINUC_READNUM_TYPES = {"scnome", "snmcseq2"}


def load_trinuc_data(config_dir, dtype_label_pairs):
    """
    Load trinuc summary TSVs and return a mapping:
        {dataset_label: DataFrame(index=CellID, columns=[noncpg, endo, exo])}

    For methods with per-read rows (R1/R2 encoded as SRR_1 / SRR_2),
    rows are grouped by prefix (strip trailing _1 or _2) and averaged.
    """
    config_dir = pathlib.Path(config_dir)
    result = {}

    for dtype, label in dtype_label_pairs:
        rel = TRINUC_FILE_MAP.get(dtype)
        if rel is None:
            continue
        path = config_dir / rel
        if not path.exists():
            print(f"  [trinuc] SKIP {label}: {path} not found")
            continue

        tdf = pd.read_csv(path, sep="\t")
        tdf = tdf.rename(columns={"sample": "CellID"})

        if dtype in TRINUC_READNUM_TYPES:
            # Strip _1 / _2 suffix and average R1+R2 per cell
            tdf["CellID"] = tdf["CellID"].str.replace(r'_\d+$', '', regex=True)
            tdf = tdf.groupby("CellID", as_index=True)[["noncpg", "endo", "exo"]].mean()
        else:
            tdf = tdf.set_index("CellID")[["noncpg", "endo", "exo"]]

        result[label] = tdf
        print(f"  [trinuc] {label:20s}  {len(tdf)} cells  ({path.name})")

    return result


def merge_trinuc_into_hdf(hdf, trinuc_data):
    """
    Add NonCG_rate and GCH_rate columns to hdf by matching Dataset + CellID.
    """
    if "NonCG_rate" not in hdf.columns:
        hdf["NonCG_rate"] = np.nan
    if "GCH_rate" not in hdf.columns:
        hdf["GCH_rate"] = np.nan

    for label, tdf in trinuc_data.items():
        mask = hdf["Dataset"] == label
        if not mask.any():
            continue
        matched = hdf.loc[mask, "CellID"].map(tdf["noncpg"])
        hdf.loc[mask, "NonCG_rate"] = matched.values

        if "exo" in tdf.columns:
            matched_exo = hdf.loc[mask, "CellID"].map(tdf["exo"])
            hdf.loc[mask, "GCH_rate"] = matched_exo.values

    return hdf


# ═══════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════

COLOR_PALETTE = [
    "#4C78A8", "#E45756", "#72B7B2", "#59A14F",
    "#F28E2B", "#B279A2", "#FF9DA7", "#9C755F", "#BAB0AC",
]


def _violin_box_plot(ax, data_series_list, labels, colors, title):
    """Draw violin + box + scatter for a list of (values, label, color) triples."""
    positions = list(range(1, len(data_series_list) + 1))

    valid = [(s, lbl, c, pos)
             for s, lbl, c, pos in zip(data_series_list, labels, colors, positions)
             if len(s) > 0]

    if not valid:
        ax.set_title(title, fontsize=10)
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="#aaaaaa")
        return

    vdata  = [s for s, *_ in valid]
    vlabels = [lbl for _, lbl, *_ in valid]
    vcolors = [c for _, _, c, *_ in valid]
    vpos   = list(range(1, len(valid) + 1))

    parts = ax.violinplot(vdata, positions=vpos,
                          showmeans=False, showmedians=False, showextrema=False)
    for body, color in zip(parts["bodies"], vcolors):
        body.set_facecolor(color)
        body.set_alpha(0.25)

    bp = ax.boxplot(vdata, positions=vpos, widths=0.18, patch_artist=True,
                    showfliers=False, zorder=3)
    for patch, color in zip(bp["boxes"], vcolors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    for element in ("whiskers", "caps", "medians"):
        for line in bp[element]:
            line.set_color("#333333")
            line.set_linewidth(1.1)

    for pos, vals, color in zip(vpos, vdata, vcolors):
        jitter = np.random.default_rng(42).normal(0, 0.04, size=len(vals))
        ax.scatter(pos + jitter, vals, c=color, alpha=0.35, s=6,
                   zorder=2, edgecolors="none")
        med = np.median(vals)
        ax.annotate(f"{med:.3g}", xy=(pos, med),
                    xytext=(pos + 0.28, med), fontsize=6.5, color="#555555",
                    arrowprops=dict(arrowstyle="-", color="#dddddd", lw=0.5))

    ax.set_title(title, fontsize=10)
    ax.set_xticks(vpos)
    ax.set_xticklabels(vlabels, fontsize=8, rotation=25, ha="right")
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_metric_group(hdf, group_metrics, group_title, pdf, dataset_order=None):
    """
    Plot a group of canonical metrics for all datasets that have that column.
    One page per group.
    """
    if dataset_order is None:
        dataset_order = hdf["Dataset"].unique().tolist()

    color_map = {lbl: COLOR_PALETTE[i % len(COLOR_PALETTE)]
                 for i, lbl in enumerate(dataset_order)}

    valid_metrics = []
    for col, display in group_metrics:
        if col in hdf.columns:
            has_data = False
            for lbl in dataset_order:
                sub = pd.to_numeric(
                    hdf.loc[hdf["Dataset"] == lbl, col], errors="coerce"
                ).dropna()
                if len(sub) > 0:
                    has_data = True
                    break
            if has_data:
                valid_metrics.append((col, display))

    if not valid_metrics:
        return

    n = len(valid_metrics)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 5.5), squeeze=False)
    axes = axes[0]
    fig.suptitle(group_title, fontsize=14, fontweight="bold", y=1.01)

    for ax, (col, display) in zip(axes, valid_metrics):
        series_list, lbls, colors = [], [], []
        for lbl in dataset_order:
            s = pd.to_numeric(
                hdf.loc[hdf["Dataset"] == lbl, col], errors="coerce"
            ).dropna().values
            if len(s) > 0:
                series_list.append(s)
                lbls.append(lbl)
                colors.append(color_map[lbl])
        _violin_box_plot(ax, series_list, lbls, colors, display)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def build_summary_table_page(hdf, all_metrics, dataset_order, pdf):
    """Append a summary statistics table page to the PDF."""
    rows = []
    for col, display in all_metrics:
        if col not in hdf.columns:
            continue
        row = {"Metric": display}
        has_any = False
        for lbl in dataset_order:
            s = pd.to_numeric(
                hdf.loc[hdf["Dataset"] == lbl, col], errors="coerce"
            ).dropna()
            if len(s) > 0:
                has_any = True
                row[f"{lbl}\nMedian"] = f"{s.median():.3g}"
                row[f"{lbl}\nMean"]   = f"{s.mean():.3g}"
            else:
                row[f"{lbl}\nMedian"] = "—"
                row[f"{lbl}\nMean"]   = "—"
        if has_any:
            rows.append(row)

    if not rows:
        return

    fig, ax = plt.subplots(figsize=(max(10, 2.5 * len(dataset_order)), 8))
    ax.axis("off")
    ax.set_title("Summary table (median / mean per method)", fontsize=12,
                 fontweight="bold", pad=16)

    col_labels = ["Metric"] + [f"{lbl}\nMedian" for lbl in dataset_order] + \
                              [f"{lbl}\nMean"   for lbl in dataset_order]
    cell_text  = [[r.get(c, "—") for c in col_labels] for r in rows]

    tbl = ax.table(cellText=cell_text, colLabels=col_labels, loc="center",
                   cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)
    tbl.scale(1.0, 1.35)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#4C78A8")
            cell.set_text_props(color="white", fontweight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f0f4f8")

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════
# Summary CSV
# ═══════════════════════════════════════════════════════════════════

def build_summary_csv(hdf, all_metrics, dataset_order):
    rows = []
    for col, display in all_metrics:
        if col not in hdf.columns:
            continue
        row = {"Metric": col, "MetricLabel": display}
        for lbl in dataset_order:
            s = pd.to_numeric(
                hdf.loc[hdf["Dataset"] == lbl, col], errors="coerce"
            ).dropna()
            row[f"{lbl}_N"]      = len(s)
            row[f"{lbl}_Median"] = s.median() if len(s) else None
            row[f"{lbl}_Mean"]   = s.mean()   if len(s) else None
            row[f"{lbl}_Std"]    = s.std()    if len(s) else None
        rows.append(row)
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Unified QC comparison across all benchmark methods"
    )
    parser.add_argument("--config", required=True,
                        help="CSV/TSV config file with columns: path,label,type")
    parser.add_argument("--output", default="all_methods_qc",
                        help="Output file prefix (default: all_methods_qc)")
    parser.add_argument("--mapq", type=int, default=30,
                        help="MAPQ label for titles (default: 30)")
    args = parser.parse_args()

    sep = "\t" if args.config.endswith(".tsv") else ","
    config = pd.read_csv(args.config, sep=sep)
    if not {"path", "label", "type"}.issubset(config.columns):
        raise ValueError("Config must have columns: path, label, type")

    print("Loading datasets:")
    loaded_rows = []
    harm_dfs = []
    for _, row in config.iterrows():
        df = load_dataset(row["path"], row["label"], row["type"])
        if df is None:
            continue
        loaded_rows.append(row)
        h = harmonize_dataset(df, row["type"], row["label"])
        harm_dfs.append(h)

    if not harm_dfs:
        print("ERROR: No datasets could be loaded. Check paths in config.")
        return

    skipped = len(config) - len(loaded_rows)
    if skipped:
        print(f"\nWARNING: {skipped} dataset(s) skipped (files not found).")

    hdf = pd.concat(harm_dfs, ignore_index=True)
    loaded_config = pd.DataFrame(loaded_rows)
    dataset_order = loaded_config["label"].tolist()

    # Classify loaded datasets by modality
    methyl_datasets = [r["label"] for _, r in loaded_config.iterrows()
                       if r["type"] in METHYL_TYPES]
    hic_datasets    = [r["label"] for _, r in loaded_config.iterrows()
                       if r["type"] in HIC_TYPES]

    # Merge pre-computed trinuc noncpg/exo rates from summary/trinuc/ folder
    config_dir = pathlib.Path(args.config).parent
    dtype_label_pairs = [(r["type"], r["label"]) for _, r in loaded_config.iterrows()]
    print("\nLoading trinuc conversion data:")
    trinuc_data = load_trinuc_data(config_dir, dtype_label_pairs)
    if trinuc_data:
        hdf = merge_trinuc_into_hdf(hdf, trinuc_data)

    print(f"\nTotal cells after harmonization: {len(hdf)}")
    for lbl in dataset_order:
        n = (hdf["Dataset"] == lbl).sum()
        print(f"  {lbl}: {n}")

    # Save harmonized CSV
    harm_path = f"{args.output}.harmonized.csv"
    hdf.to_csv(harm_path, index=False)
    print(f"\nSaved harmonized data: {harm_path}")

    all_metrics = BASE_ALIGNMENT + METHYLATION + HIC_CONTACTS

    # Summary CSV
    summary = build_summary_csv(hdf, all_metrics, dataset_order)
    summary_path = f"{args.output}.summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"Saved summary: {summary_path}")

    print("\n--- Median overview (MapQ30Rate / DupRate / mCG_Rate) ---")
    for lbl in dataset_order:
        sub = hdf[hdf["Dataset"] == lbl]
        def _med(col):
            return pd.to_numeric(sub[col], errors="coerce").median() if col in sub.columns else float("nan")
        mq  = _med("MapQ30Rate")
        dup = _med("DupRate")
        mcg = _med("mCG_Rate")
        mq_str  = f"{mq:5.1f}%" if pd.notna(mq)  else "   N/A"
        dup_str = f"{dup:5.1f}%" if pd.notna(dup) else "   N/A"
        mcg_str = f"{mcg:.1f}%"  if pd.notna(mcg) else "N/A"
        print(f"  {lbl:22s}  MapQ30Rate={mq_str}  DupRate={dup_str}  mCG={mcg_str}")

    # ── PDF 1: Base alignment (all methods) ──────────────────────────
    base_pdf = f"{args.output}.base_alignment.pdf"
    with PdfPages(base_pdf) as pdf:
        plot_metric_group(hdf, BASE_ALIGNMENT,
                          "Base Alignment QC — All Methods",
                          pdf, dataset_order)
        build_summary_table_page(hdf, BASE_ALIGNMENT, dataset_order, pdf)
    print(f"\nSaved: {base_pdf}")

    # ── PDF 2: Methylation (bisulfite methods only) ───────────────────
    meth_hdf = hdf[hdf["Dataset"].isin(methyl_datasets)]
    meth_pdf = f"{args.output}.methylation.pdf"
    with PdfPages(meth_pdf) as pdf:
        plot_metric_group(meth_hdf, METHYLATION,
                          "Methylation QC — Bisulfite Methods",
                          pdf, methyl_datasets)
        build_summary_table_page(meth_hdf, METHYLATION, methyl_datasets, pdf)
    print(f"Saved: {meth_pdf}")

    # ── PDF 3: HiC contacts (HiC methods only) ───────────────────────
    hic_hdf = hdf[hdf["Dataset"].isin(hic_datasets)]
    hic_pdf = f"{args.output}.hic_contacts.pdf"
    with PdfPages(hic_pdf) as pdf:
        plot_metric_group(hic_hdf, HIC_CONTACTS,
                          "Hi-C Contact QC — 3C Methods",
                          pdf, hic_datasets)
        build_summary_table_page(hic_hdf, HIC_CONTACTS, hic_datasets, pdf)
    print(f"Saved: {hic_pdf}")

    print("\nDone.")


if __name__ == "__main__":
    main()
