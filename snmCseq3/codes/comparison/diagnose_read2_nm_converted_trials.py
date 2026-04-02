#!/usr/bin/env python3
"""Diagnose read2 ``NM:i`` vs converted-genome recomputation (no MD for heuristics).

1. **Classification** of each read2: is ``NM`` among the four CT/GA×orientation trials? Does
   ``min(trials)`` equal ``NM``? Does ``MD``-based recompute agree with ``NM``?

2. **Base-level samples** for reads where ``NM`` is not in the trial set (or ``NM`` ≠ ``MD``):
   print aligned columns with SAM base, ``MD`` reference base, CT-genome base, GA-genome base,
   and converted-query bases (both orientations).

3. **Aggregation tests** (still no ``MD`` for picking among trials): compare accuracy vs ``NM`` for
   ``min(all 4)``, ``min(CT only)``, ``min(GA only)``, and **unique minimum** (if the global minimum
   distance appears in exactly one trial, use it; otherwise fall back to ``min(all)`` — a
   tie-breaker variant).

Example::

  python diagnose_read2_nm_converted_trials.py \\
    /gpfs/.../mm10_bismark/Bisulfite_Genome \\
    /gpfs/.../mm10/mm10.fa \\
    alignments.bam --max-reads 8000 --detail-n 15
"""

from __future__ import annotations

import argparse
import collections
import os
import sys
from typing import Iterable

import numpy as np
import pysam

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bisulfite_corrected_mismatch import (  # noqa: E402
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
    count_nm_style_edit_distance_from_md,
    _pbat_converted_query_variants,
)


def _labeled_trials(
    read, fa_ct: pysam.FastaFile, fa_ga: pysam.FastaFile
) -> list[tuple[str, int]]:
    c_ct = bisulfite_converted_contig_name(fa_ct, read.reference_name, "CT")
    c_ga = bisulfite_converted_contig_name(fa_ga, read.reference_name, "GA")
    if c_ct is None or c_ga is None:
        return []
    out: list[tuple[str, int]] = []
    for qtag, qconv in _pbat_converted_query_variants(read):
        d_ct = count_nm_style_edit_distance_converted_explicit(
            read, fa_ct, ref_contig=c_ct, query_converted=qconv
        )
        d_ga = count_nm_style_edit_distance_converted_explicit(
            read, fa_ga, ref_contig=c_ga, query_converted=qconv
        )
        if d_ct >= 0:
            out.append((f"CT_{qtag}", d_ct))
        if d_ga >= 0:
            out.append((f"GA_{qtag}", d_ga))
    return out


def _min_for_labels(
    trials: list[tuple[str, int]], prefixes: Iterable[str]
) -> int | None:
    wanted = [d for lab, d in trials if any(lab.startswith(p) for p in prefixes)]
    return min(wanted) if wanted else None


def _second_smallest_dist(trials: list[tuple[str, int]]) -> int | None:
    ds = sorted({d for _, d in trials})
    if not ds:
        return None
    return ds[1] if len(ds) > 1 else ds[0]


def _detail_rows(
    read,
    fa_ct: pysam.FastaFile,
    fa_ga: pysam.FastaFile,
    fa_genome: pysam.FastaFile | None,
    limit: int = 12,
) -> list[str]:
    chrom = read.reference_name
    c_ct = bisulfite_converted_contig_name(fa_ct, chrom, "CT")
    c_ga = bisulfite_converted_contig_name(fa_ga, chrom, "GA")
    if c_ct is None or c_ga is None:
        return []
    qs = read.query_sequence
    lines: list[str] = []
    qv = list(_pbat_converted_query_variants(read))
    n = 0
    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except (ValueError, AttributeError):
        return ["(no aligned pairs with_seq)"]
    for t in pairs:
        if len(t) != 3:
            continue
        q, rpos, ref_md = t
        if q is None or rpos is None or ref_md is None:
            continue
        sam_b = qs[q].upper()
        rmd = ref_md.upper()
        ct_b = fa_ct.fetch(c_ct, rpos, rpos + 1).upper()
        ga_b = fa_ga.fetch(c_ga, rpos, rpos + 1).upper()
        g_b = ""
        if fa_genome is not None:
            try:
                g_b = fa_genome.fetch(chrom, rpos, rpos + 1).upper()
            except (ValueError, KeyError):
                g_b = "?"
        q_ct_seq = qv[0][1][q].upper()
        q_ct_rc = qv[1][1][q].upper()
        mismatch_md = sam_b != rmd
        mismatch_ct_seq = q_ct_seq != ct_b
        lines.append(
            f"  q={q} rpos={rpos} SAM={sam_b} MDref={rmd} mm10={g_b or '-'} "
            f"CTref={ct_b} GAref={ga_b} qCTseq={q_ct_seq} qCTrc={q_ct_rc} "
            f"diffMD={mismatch_md} diffCTseq={mismatch_ct_seq}"
        )
        n += 1
        if n >= limit:
            break
    return lines if lines else ["(no match columns)"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bisulfite_genome_dir", help="Bisulfite_Genome with CT_conversion/ GA_conversion/")
    ap.add_argument("genomic_fasta", help="Unconverted mm10.fa (optional detail); use '' to skip")
    ap.add_argument("bam")
    ap.add_argument("--max-reads", type=int, default=8000)
    ap.add_argument("--detail-n", type=int, default=12, help="reads to dump base columns for")
    args = ap.parse_args()

    base = args.bisulfite_genome_dir.rstrip("/")
    ct_fa = f"{base}/CT_conversion/genome_mfa.CT_conversion.fa"
    ga_fa = f"{base}/GA_conversion/genome_mfa.GA_conversion.fa"
    for p in (ct_fa, ga_fa):
        if not os.path.isfile(p):
            print(f"ERROR: missing {p}", file=sys.stderr)
            sys.exit(1)

    fa_genome: pysam.FastaFile | None = None
    if args.genomic_fasta and args.genomic_fasta.strip():
        if not os.path.isfile(args.genomic_fasta):
            print(f"WARN: genomic fasta missing, skip mm10 column: {args.genomic_fasta}", file=sys.stderr)
        else:
            fa_genome = pysam.FastaFile(args.genomic_fasta)

    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Counters
    n = 0
    nm_in_trials = 0
    min_eq_nm = 0
    min_ne_nm_but_nm_in = 0
    nm_not_in_trials = 0
    nm_ne_md = 0
    nm_eq_md = 0

    strat_min_all = strat_min_ct = strat_min_ga = strat_second = 0

    # min(all) mismatch reasons
    cause_min_too_low = 0  # min < NM
    cause_min_too_high = 0  # min > NM

    detail_budget = args.detail_n
    printed_nm_not_in = 0
    printed_nm_ne_md = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.is_read2 or not read.has_tag("NM"):
            continue

        trials = _labeled_trials(read, fa_ct, fa_ga)
        if len(trials) != 4:
            continue

        nm = int(read.get_tag("NM"))
        mdv = count_nm_style_edit_distance_from_md(read)
        dists = [d for _, d in trials]
        m_all = min(dists)
        m_ct = _min_for_labels(trials, ("CT_",))
        m_ga = _min_for_labels(trials, ("GA_",))
        m_2nd = _second_smallest_dist(trials)

        trial_vals = set(dists)
        if nm in trial_vals:
            nm_in_trials += 1
        else:
            nm_not_in_trials += 1

        if m_all == nm:
            min_eq_nm += 1
        elif nm in trial_vals:
            min_ne_nm_but_nm_in += 1
            if m_all < nm:
                cause_min_too_low += 1
            else:
                cause_min_too_high += 1

        if m_all == nm:
            strat_min_all += 1
        if m_ct == nm:
            strat_min_ct += 1
        if m_ga == nm:
            strat_min_ga += 1
        if m_2nd == nm:
            strat_second += 1

        if mdv >= 0:
            if mdv != nm:
                nm_ne_md += 1
            else:
                nm_eq_md += 1

        if nm not in trial_vals and printed_nm_not_in < detail_budget:
            print("\n--- read2 NM not in {CT,GA}×2 trials ---")
            print(
                f"qname={read.query_name} chrom={read.reference_name} pos={read.reference_start} "
                f"rev={read.is_reverse} NM={nm} MDrec={mdv} min4={m_all} trials={trials}"
            )
            for row in _detail_rows(read, fa_ct, fa_ga, fa_genome, limit=15):
                print(row)
            printed_nm_not_in += 1

        if mdv >= 0 and mdv != nm and printed_nm_ne_md < min(5, detail_budget):
            print("\n--- read2 NM != MD recompute (tag inconsistency) ---")
            print(
                f"qname={read.query_name} NM={nm} MDrec={mdv} min4={m_all} trials={trials}"
            )
            for row in _detail_rows(read, fa_ct, fa_ga, fa_genome, limit=10):
                print(row)
            printed_nm_ne_md += 1

        n += 1
        if n >= args.max_reads:
            break

    bam.close()
    fa_ct.close()
    fa_ga.close()
    if fa_genome:
        fa_genome.close()

    print("=== Read2 diagnosis (primary, mapped) ===")
    print(f"reads_used\t{n}")
    print(f"NM appears in at least one of 4 trials\t{nm_in_trials}\t({nm_in_trials/n:.4f})")
    print(f"min(all4) == NM\t{min_eq_nm}\t({min_eq_nm/n:.4f})")
    print(
        f"min(all4) != NM but NM in trial set\t{min_ne_nm_but_nm_in}\t"
        f"({min_ne_nm_but_nm_in/n:.4f})  [min<NM: {cause_min_too_low}, min>NM: {cause_min_too_high}]"
    )
    print(f"NM not in any trial\t{nm_not_in_trials}\t({nm_not_in_trials/n:.4f})")
    print()
    print("NM vs MD (informational):")
    print(f"  NM == MDrec\t{nm_eq_md}\tNM != MDrec\t{nm_ne_md}")
    print()
    print("Heuristic vs NM (fraction rec == NM):")
    print(f"  min(all 4 trials)\t{strat_min_all/n:.4f}")
    print(f"  min(CT trials only)\t{strat_min_ct/n:.4f}")
    print(f"  min(GA trials only)\t{strat_min_ga/n:.4f}")
    print(
        f"  second-smallest distinct distance among 4\t{strat_second/n:.4f}  "
        "(using fewer trials as 'not the minimum' — usually worse)"
    )
    print()
    print(
        "Interpretation: If NM is in the trial set but min(all4) != NM, the error is "
        "**wrong index/orientation choice** (another trial scores lower spuriously). "
        "If NM is not in any trial, our four converted strings do not match BWA’s query "
        "encoding and/or NM/MD disagree."
    )


if __name__ == "__main__":
    main()
