#!/usr/bin/env python3
"""
NM vs **best-effort recomputed** edit distance (read1 / read2), plus MAPQ stratification vs yap.

**Bhmem:** ``recompute_nm_from_converted_genomes_pbat(..., use_md_fallback=True)`` — maximizes
agreement with ``NM:i`` in benchmarks (CT/GA trials + ``MD`` when needed).

**Yap/Bismark:** ``recompute_nm_style_from_md`` when ``MD:Z`` is present (typically matches ``NM``).

**Cross-pipeline comparable** mismatch: ``count_cross_pipeline_comparable_edit_distance`` on **mm10.fa**
(same call on both BAMs) for “which alignment is better under one definition”.

Strata (from ``joined.tsv`` ``mapq_bhmem``, ``mapq_yap``):

- **all:** subsampled rows for overview scatter plots.
- **yap_high_bhmem_low:** ``mapq_bhmem < 30`` and ``mapq_yap > 30`` (yap confident, bhmem not).
- **diff_locus:** same MAPQ filter plus **not** same locus (different ``RNAME`` or ``|start| > 50``).

Outputs: ``*.per_read.tsv``, ``*.summary.txt``, ``*.plots.pdf``.

Example::

  python report_nm_replicate_mapq_stratified.py \\
    mapq_comparison/SRR21549292/mapq_comparison.joined.tsv \\
    04.bhmem_bam/SRR21549292.bhmem.bam \\
    /path/to/SRR21549292.3C.sorted.bam \\
    /gpfs/.../mm10/mm10.fa \\
    /gpfs/.../mm10_bismark/Bisulfite_Genome \\
    -o mapq_comparison/SRR21549292/nm_replicate_mapq_stratified \\
    --max-rows-all 25000 --max-rows-stratum 40000
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    plt = None
    PdfPages = None

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bisulfite_corrected_mismatch import (  # noqa: E402
    count_cross_pipeline_comparable_edit_distance,
    recompute_nm_from_converted_genomes_pbat,
    recompute_nm_style_from_md,
)


def parse_yap_qname(qname: str):
    parts = qname.split("_")
    if len(parts) < 2:
        return None, None
    base, strand = parts[0], parts[1]
    if strand.startswith("1"):
        return base, True
    if strand.startswith("2"):
        return base, False
    return None, None


def collect_bhmem_primary(bam_path: str, keys: set) -> dict:
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not (read.flag & 64) and not (read.flag & 128):
                continue
            is_r1 = bool(read.flag & 64)
            key = (read.query_name, is_r1)
            if key in keys:
                out[key] = read
    return out


def collect_yap_best(bam_path: str, keys: set) -> dict:
    by_key = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            base, is_r1 = parse_yap_qname(read.query_name)
            if base is None:
                continue
            key = (base, is_r1)
            if key not in keys:
                continue
            parts = read.query_name.split("_")
            is_prim = "-" not in parts[-1]
            mq = read.mapping_quality
            by_key[key].append((read, mq, is_prim))
    best = {}
    for key, cands in by_key.items():
        prim = [x for x in cands if x[2]]
        pool = prim if prim else cands
        best[key] = max(pool, key=lambda x: x[1])[0]
    return best


def load_joined(path: str) -> list[dict]:
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(
                {
                    "base_id": row["base_id"],
                    "is_r1": int(row["is_r1"]),
                    "mapq_bhmem": int(row["mapq_bhmem"]),
                    "mapq_yap": int(row["mapq_yap"]),
                }
            )
    return rows


def nm_replicate_bhmem(read, fa_ct, fa_ga) -> int:
    return recompute_nm_from_converted_genomes_pbat(
        read, fa_ct, fa_ga, use_md_fallback=True
    )


def nm_replicate_yap(read) -> int:
    if not read.has_tag("MD"):
        return -1
    return recompute_nm_style_from_md(read)


def same_locus(br, yr, max_delta: int) -> bool:
    if br is None or yr is None:
        return False
    if br.reference_name != yr.reference_name:
        return False
    return abs(br.reference_start - yr.reference_start) <= max_delta


def _scatter(ax, nm, rec, title, xlab="NM tag", ylab="Recomputed edit"):
    nm = np.asarray(nm, dtype=np.float64)
    rec = np.asarray(rec, dtype=np.float64)
    m = (nm >= 0) & (rec >= 0)
    if not np.any(m):
        ax.set_title(title + "\n(no data)")
        return
    ax.scatter(nm[m], rec[m], alpha=0.15, s=6, c="#1a5276", edgecolors="none", rasterized=True)
    hi = max(float(nm[m].max()), float(rec[m].max()))
    ax.plot([0, hi], [0, hi], "k--", lw=0.8, alpha=0.4)
    ax.set_xlim(left=-0.5)
    ax.set_ylim(bottom=-0.5)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if np.sum(m) >= 2:
        r = np.corrcoef(nm[m], rec[m])[0, 1]
        frac = float(np.mean(nm[m] == rec[m]))
        ax.set_title(f"{title}\nn={int(np.sum(m))}  frac(eq)={frac:.3f}  r={r:.3f}")
    else:
        ax.set_title(title)


def _hist_pair(ax, a, b, la, lb, title, bins=35):
    a = np.asarray([x for x in a if x >= 0], dtype=np.int64)
    b = np.asarray([x for x in b if x >= 0], dtype=np.int64)
    if len(a) == 0 and len(b) == 0:
        ax.set_title(title + "\n(no data)")
        return
    hi = 1
    if len(a):
        hi = max(hi, int(a.max()))
    if len(b):
        hi = max(hi, int(b.max()))
    edges = np.arange(-0.5, min(hi, 80) + 1.5, max(1, (min(hi, 80) + 1) // bins))
    ax.hist(
        a,
        bins=edges,
        alpha=0.55,
        label=f"{la} (n={len(a)})",
        color="#1a5276",
        density=True,
    )
    ax.hist(
        b,
        bins=edges,
        alpha=0.55,
        label=f"{lb} (n={len(b)})",
        color="#c0392b",
        density=True,
    )
    ax.set_xlabel("Mismatch / edit distance")
    ax.set_ylabel("Density")
    ax.legend(fontsize=7)
    ax.set_title(title)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("joined_tsv")
    ap.add_argument("bhmem_bam")
    ap.add_argument("yap_bam")
    ap.add_argument("reference_fasta", help="Genomic mm10.fa (unified metric)")
    ap.add_argument(
        "bisulfite_genome_dir",
        help="Bisulfite_Genome with CT_conversion/ GA_conversion/ (bhmem replicate)",
    )
    ap.add_argument("-o", "--output-prefix", default="nm_replicate_mapq_stratified")
    ap.add_argument(
        "--max-rows-all",
        type=int,
        default=25_000,
        help="Max joined rows for full scatter (0 = use all; slow)",
    )
    ap.add_argument(
        "--max-rows-stratum",
        type=int,
        default=40_000,
        help="Cap stratum rows after filter (0 = no cap)",
    )
    ap.add_argument(
        "--bhmem-mapq-lt",
        type=int,
        default=30,
        help="Keep rows with bhmem MAPQ strictly less than this (default 30)",
    )
    ap.add_argument(
        "--yap-mapq-gt",
        type=int,
        default=30,
        help="Keep rows with yap MAPQ strictly greater than this (default 30)",
    )
    ap.add_argument(
        "--locus-delta",
        type=int,
        default=50,
        help="|start_bh - start_yap| <= this and same chr => same locus",
    )
    args = ap.parse_args()

    base_bg = args.bisulfite_genome_dir.rstrip("/")
    ct_fa = f"{base_bg}/CT_conversion/genome_mfa.CT_conversion.fa"
    ga_fa = f"{base_bg}/GA_conversion/genome_mfa.GA_conversion.fa"
    for p in (ct_fa, ga_fa):
        if not os.path.isfile(p):
            print(f"ERROR: missing {p}", file=sys.stderr)
            sys.exit(1)

    rows = load_joined(args.joined_tsv)
    print(f"Loaded {len(rows)} joined rows", file=sys.stderr)

    # keys for all-sample plot
    n_all = args.max_rows_all if args.max_rows_all else len(rows)
    rows_all = rows[: min(n_all, len(rows))]
    keys_all = {(r["base_id"], bool(r["is_r1"])) for r in rows_all}

    stratum_rows = [
        r
        for r in rows
        if r["mapq_bhmem"] < args.bhmem_mapq_lt and r["mapq_yap"] > args.yap_mapq_gt
    ]
    if args.max_rows_stratum and len(stratum_rows) > args.max_rows_stratum:
        stratum_rows = stratum_rows[: args.max_rows_stratum]

    keys_stratum = {(r["base_id"], bool(r["is_r1"])) for r in stratum_rows}
    keys_union = keys_all | keys_stratum

    print(f"Collecting BAM records for {len(keys_union)} keys...", file=sys.stderr)
    bh = collect_bhmem_primary(args.bhmem_bam, keys_union)
    yp = collect_yap_best(args.yap_bam, keys_union)
    print(f"  bhmem {len(bh)}  yap {len(yp)}", file=sys.stderr)

    fa = pysam.FastaFile(args.reference_fasta)
    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)

    per_read = f"{args.output_prefix}.per_read.tsv"
    summary_path = f"{args.output_prefix}.summary.txt"

    def process_row(rec: dict) -> dict | None:
        key = (rec["base_id"], bool(rec["is_r1"]))
        br, yr = bh.get(key), yp.get(key)
        if br is None or yr is None:
            return None
        nb = int(br.get_tag("NM")) if br.has_tag("NM") else -1
        ny = int(yr.get_tag("NM")) if yr.has_tag("NM") else -1
        rb = nm_replicate_bhmem(br, fa_ct, fa_ga)
        ry = nm_replicate_yap(yr)
        ub = count_cross_pipeline_comparable_edit_distance(br, fa)
        uy = count_cross_pipeline_comparable_edit_distance(yr, fa)
        sl = same_locus(br, yr, args.locus_delta)
        return {
            "base_id": rec["base_id"],
            "is_r1": rec["is_r1"],
            "mapq_bhmem": rec["mapq_bhmem"],
            "mapq_yap": rec["mapq_yap"],
            "nm_bhmem": nb,
            "nm_yap": ny,
            "replicate_bhmem": rb,
            "replicate_yap": ry,
            "unified_bhmem": ub,
            "unified_yap": uy,
            "chr_bhmem": br.reference_name,
            "pos_bhmem": br.reference_start,
            "chr_yap": yr.reference_name,
            "pos_yap": yr.reference_start,
            "same_locus": int(sl),
        }

    written = []
    with open(per_read, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(
            [
                "base_id",
                "is_r1",
                "mapq_bhmem",
                "mapq_yap",
                "nm_bhmem",
                "nm_yap",
                "replicate_bhmem",
                "replicate_yap",
                "unified_bhmem",
                "unified_yap",
                "chr_bhmem",
                "pos_bhmem",
                "chr_yap",
                "pos_yap",
                "same_locus",
            ]
        )
        for rec in rows:
            if (rec["base_id"], bool(rec["is_r1"])) not in keys_union:
                continue
            o = process_row(rec)
            if o is None:
                continue
            w.writerow(
                [
                    o["base_id"],
                    o["is_r1"],
                    o["mapq_bhmem"],
                    o["mapq_yap"],
                    o["nm_bhmem"],
                    o["nm_yap"],
                    o["replicate_bhmem"],
                    o["replicate_yap"],
                    o["unified_bhmem"],
                    o["unified_yap"],
                    o["chr_bhmem"],
                    o["pos_bhmem"],
                    o["chr_yap"],
                    o["pos_yap"],
                    o["same_locus"],
                ]
            )
            written.append(o)

    fa.close()
    fa_ct.close()
    fa_ga.close()

    by_key = {(o["base_id"], o["is_r1"]): o for o in written}

    def pick_fast(keys_set, pred):
        out = []
        for rec in rows:
            key = (rec["base_id"], bool(rec["is_r1"]))
            if key not in keys_set:
                continue
            o = by_key.get(key)
            if o is None:
                continue
            if pred(o):
                out.append(o)
        return out

    keys_all_f = {(r["base_id"], bool(r["is_r1"])) for r in rows_all}
    keys_stratum_f = {(r["base_id"], bool(r["is_r1"])) for r in stratum_rows}

    def split_r12(seq):
        r1 = [x for x in seq if x["is_r1"]]
        r2 = [x for x in seq if not x["is_r1"]]
        return r1, r2

    all_o = [by_key[k] for k in keys_all_f if k in by_key]
    strat_o = pick_fast(keys_stratum_f, lambda o: True)
    diff_o = [o for o in strat_o if not o["same_locus"]]

    lines = []
    lines.append(f"n_joined_file\t{len(rows)}")
    lines.append(f"n_written_per_read\t{len(written)}")
    lines.append(f"n_all_plot_keys\t{len(all_o)}")
    lines.append(f"n_stratum_mapq\t{len(strat_o)}")
    lines.append(f"n_stratum_diff_locus\t{len(diff_o)}")
    lines.append("")

    def stat_block(name, seq):
        if not seq:
            lines.append(f"## {name}\n(empty)\n")
            return
        pu = [
            (o["unified_bhmem"], o["unified_yap"])
            for o in seq
            if o["unified_bhmem"] >= 0 and o["unified_yap"] >= 0
        ]
        pr = [
            (o["replicate_bhmem"], o["replicate_yap"])
            for o in seq
            if o["replicate_bhmem"] >= 0 and o["replicate_yap"] >= 0
        ]
        lines.append(f"## {name}")
        lines.append(f"n\t{len(seq)}")
        if pu:
            ub = np.array([a for a, _ in pu], dtype=np.float64)
            uy = np.array([b for _, b in pu], dtype=np.float64)
            lines.append(
                f"mean_unified_yap_minus_bhmem\t{float(np.mean(uy - ub)):.6f}"
            )
            lines.append(
                f"median_unified_yap_minus_bhmem\t{float(np.median(uy - ub)):.6f}"
            )
            lines.append(
                f"fraction_unified_yap_lt_bhmem\t{float(np.mean(uy < ub)):.6f}"
            )
        if pr:
            rb = np.array([a for a, _ in pr], dtype=np.float64)
            ry = np.array([b for _, b in pr], dtype=np.float64)
            lines.append(
                f"mean_replicate_yap_minus_bhmem\t{float(np.mean(ry - rb)):.6f}"
            )
            lines.append(
                f"fraction_replicate_yap_lt_bhmem\t{float(np.mean(ry < rb)):.6f}"
            )
        lines.append("")

    stat_block("stratum_mapq_bhmem_lt_30_yap_gt_30", strat_o)
    stat_block("stratum_plus_diff_locus", diff_o)

    with open(summary_path, "w") as sf:
        sf.write("\n".join(lines))

    print(f"Wrote {per_read}", file=sys.stderr)
    print(f"Wrote {summary_path}", file=sys.stderr)

    if plt is None:
        print("WARNING: no matplotlib, skip PDF", file=sys.stderr)
        return

    pdf_path = f"{args.output_prefix}.plots.pdf"
    with PdfPages(pdf_path) as pdf:
        # Page 1: NM vs replicate, R1/R2, bhmem + yap
        fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.5))
        a1, a2 = split_r12(all_o)
        _scatter(
            axes[0, 0],
            [x["nm_bhmem"] for x in a1],
            [x["replicate_bhmem"] for x in a1],
            "Bhmem read1",
        )
        _scatter(
            axes[0, 1],
            [x["nm_bhmem"] for x in a2],
            [x["replicate_bhmem"] for x in a2],
            "Bhmem read2",
        )
        _scatter(
            axes[1, 0],
            [x["nm_yap"] for x in a1],
            [x["replicate_yap"] for x in a1],
            "Yap read1",
        )
        _scatter(
            axes[1, 1],
            [x["nm_yap"] for x in a2],
            [x["replicate_yap"] for x in a2],
            "Yap read2",
        )
        fig.suptitle(
            "NM tag vs recomputed edit (bhmem: conv+MD; yap: MD)\n"
            f"subsample n≈{len(all_o)}",
            fontsize=10,
        )
        fig.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 2: stratum — unified and replicate histograms
        fig, axes = plt.subplots(2, 2, figsize=(9.5, 8))
        rb = [o["replicate_bhmem"] for o in strat_o]
        ry = [o["replicate_yap"] for o in strat_o]
        ub = [o["unified_bhmem"] for o in strat_o]
        uy = [o["unified_yap"] for o in strat_o]
        _hist_pair(
            axes[0, 0],
            rb,
            ry,
            "bhmem repl.",
            "yap repl.",
            "MAPQ stratum: replicates (aligner-native space)",
        )
        _hist_pair(
            axes[0, 1],
            ub,
            uy,
            "bhmem unified",
            "yap unified",
            "MAPQ stratum: unified mm10 edit (same definition)",
        )
        su = np.array(uy, dtype=np.float64) - np.array(ub, dtype=np.float64)
        axes[1, 0].hist(
            su,
            bins=min(50, max(10, int(np.ptp(su)) + 2)),
            color="#2471a3",
            edgecolor="black",
            alpha=0.85,
        )
        axes[1, 0].set_xlabel("Unified edit (yap − bhmem)")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title(
            f"Stratum n={len(strat_o)}\nnegative ⇒ bhmem higher (worse) unified edit"
        )
        axes[1, 1].axis("off")
        axes[1, 1].text(
            0.05,
            0.85,
            "Stratum: bhmem MAPQ < 30 AND yap MAPQ > 30\n\n"
            "Unified: count_cross_pipeline_comparable_edit_distance\n"
            "(XR on yap; bhmem mate-2 rule).\n\n"
            "Replicate: best NM match per pipeline.",
            fontsize=9,
            verticalalignment="top",
            fontfamily="monospace",
            transform=axes[1, 1].transAxes,
        )
        fig.suptitle("Yap-high / bhmem-low MAPQ reads", fontsize=11)
        fig.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

        # Page 3: diff locus only (if enough)
        if len(diff_o) >= 20:
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            rb = [o["replicate_bhmem"] for o in diff_o]
            ry = [o["replicate_yap"] for o in diff_o]
            ub = [o["unified_bhmem"] for o in diff_o]
            uy = [o["unified_yap"] for o in diff_o]
            _hist_pair(
                axes[0],
                ub,
                uy,
                "bhmem",
                "yap",
                f"Diff locus (n={len(diff_o)}): unified edit",
            )
            su = np.array(uy, dtype=np.float64) - np.array(ub, dtype=np.float64)
            axes[1].hist(
                su,
                bins=min(45, max(10, int(np.ptp(su)) + 2)),
                color="#884ea0",
                edgecolor="black",
                alpha=0.85,
            )
            axes[1].set_xlabel("Unified (yap − bhmem)")
            axes[1].set_ylabel("Count")
            axes[1].set_title("Difference (diff locus)")
            fig.suptitle("Subset: different locus (>50bp or different chr)", fontsize=10)
            fig.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

        # Page 4: scatter unified yap vs bhmem stratum, by mate
        if len(strat_o) >= 5:
            fig, ax = plt.subplots(figsize=(6.5, 6))
            s1 = [o for o in strat_o if o["is_r1"]]
            s2 = [o for o in strat_o if not o["is_r1"]]
            if s1:
                ax.scatter(
                    [o["unified_bhmem"] for o in s1],
                    [o["unified_yap"] for o in s1],
                    alpha=0.25,
                    s=10,
                    c="#1a5276",
                    label=f"R1 (n={len(s1)})",
                    rasterized=True,
                )
            if s2:
                ax.scatter(
                    [o["unified_bhmem"] for o in s2],
                    [o["unified_yap"] for o in s2],
                    alpha=0.25,
                    s=10,
                    c="#c0392b",
                    label=f"R2 (n={len(s2)})",
                    rasterized=True,
                )
            hi = 1
            if strat_o:
                hi = max(
                    max(o["unified_bhmem"] for o in strat_o),
                    max(o["unified_yap"] for o in strat_o),
                    1,
                )
            ax.plot([0, hi], [0, hi], "k--", lw=0.8, alpha=0.4)
            ax.set_xlabel("Unified edit bhmem")
            ax.set_ylabel("Unified edit yap")
            ax.legend()
            ax.set_title("MAPQ stratum: same cross-pipeline metric on both alignments")
            fig.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

    print(f"Wrote {pdf_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
