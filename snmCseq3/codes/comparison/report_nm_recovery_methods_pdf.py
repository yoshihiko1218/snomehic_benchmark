#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pysam
from scipy.stats import spearmanr

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("ERROR: pip install matplotlib", file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bisulfite_corrected_mismatch import (
    recompute_nm_style_from_md,
    recompute_nm_from_converted_genomes_pbat,
    recompute_nm_from_converted_genomes_pbat_no_md,
    pbat_converted_genome_trial_distances,
    count_non_md_cross_pipeline_pbat_nd_edit_distance,
)


def stats(nm, rec):
    nm = np.asarray(nm, dtype=np.int64)
    rec = np.asarray(rec, dtype=np.int64)
    m = rec >= 0
    nm = nm[m]
    rec = rec[m]
    if len(nm) == 0:
        return 0, float('nan'), float('nan'), float('nan')
    eq = float(np.mean(nm == rec))
    mad = float(np.mean(np.abs(nm - rec)))
    sp = float(spearmanr(nm, rec)[0]) if len(nm) >= 2 else float('nan')
    return len(nm), eq, mad, sp


def main():
    ap = argparse.ArgumentParser(description='Method summary table PDF for NM recovery (R1/R2).')
    ap.add_argument('bam')
    ap.add_argument('--bisulfite-genome', default='/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome')
    ap.add_argument('--mm10-fasta', default='/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa')
    ap.add_argument('--max-per-mate', type=int, default=20000)
    ap.add_argument('-o', '--output-pdf', default='nm_recovery_methods_summary.pdf')
    args = ap.parse_args()

    ct = f"{args.bisulfite_genome.rstrip('/')}/CT_conversion/genome_mfa.CT_conversion.fa"
    ga = f"{args.bisulfite_genome.rstrip('/')}/GA_conversion/genome_mfa.GA_conversion.fa"
    fa_ct = pysam.FastaFile(ct)
    fa_ga = pysam.FastaFile(ga)
    fa_mm10 = pysam.FastaFile(args.mm10_fasta)

    methods = [
        ('md_only', lambda r: recompute_nm_style_from_md(r)),
        ('conv_md_fallback', lambda r: recompute_nm_from_converted_genomes_pbat(r, fa_ct, fa_ga, use_md_fallback=True)),
        ('conv_no_md_min4', lambda r: recompute_nm_from_converted_genomes_pbat_no_md(r, fa_ct, fa_ga)),
        ('same_no_md_metric', lambda r: count_non_md_cross_pipeline_pbat_nd_edit_distance(r, fa_mm10)),
    ]

    nm = {'r1': [], 'r2': []}
    rec = {m[0]: {'r1': [], 'r2': []} for m in methods}
    trial_has_nm = {'r1': [], 'r2': []}

    bam = pysam.AlignmentFile(args.bam, 'rb')
    n1 = n2 = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.has_tag('NM'):
            continue
        mate = 'r1' if read.is_read1 else 'r2'
        if mate == 'r1' and n1 >= args.max_per_mate:
            continue
        if mate == 'r2' and n2 >= args.max_per_mate:
            continue
        t = int(read.get_tag('NM'))
        nm[mate].append(t)
        trials = pbat_converted_genome_trial_distances(read, fa_ct, fa_ga)
        trial_has_nm[mate].append(1 if t in set(trials) else 0)
        for name, fn in methods:
            rec[name][mate].append(fn(read))
        if mate == 'r1':
            n1 += 1
        else:
            n2 += 1
        if n1 >= args.max_per_mate and n2 >= args.max_per_mate:
            break
    bam.close(); fa_ct.close(); fa_ga.close(); fa_mm10.close()

    rows = []
    for mate in ('r1', 'r2'):
        th = float(np.mean(trial_has_nm[mate])) if trial_has_nm[mate] else float('nan')
        for name, _ in methods:
            n, eq, mad, sp = stats(nm[mate], rec[name][mate])
            rows.append([mate, name, n, f"{eq:.4f}", f"{mad:.4f}", f"{sp:.4f}", f"{th:.4f}"])

    # PDF table + bar chart
    with PdfPages(args.output_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(12, 5.2))
        ax.axis('off')
        col_labels = ['mate', 'method', 'n', 'frac_eq', 'mean_abs_diff', 'Spearman', 'trial_has_nm_frac']
        table = ax.table(cellText=rows, colLabels=col_labels, loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.3)
        ax.set_title('NM recovery method summary (Bhmem BAM)', fontsize=12)
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for i, mate in enumerate(('r1', 'r2')):
            labels=[]; vals=[]
            for name,_ in methods:
                labels.append(name)
                vals.append(float(stats(nm[mate], rec[name][mate])[3]))
            axes[i].bar(labels, vals, color=['#1f77b4','#2ca02c','#ff7f0e','#9467bd'])
            axes[i].set_ylim(0,1)
            axes[i].tick_params(axis='x', rotation=25)
            axes[i].set_title(f"{mate} Spearman(NM, rec)")
        fig.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

    print(args.output_pdf)


if __name__ == '__main__':
    main()
