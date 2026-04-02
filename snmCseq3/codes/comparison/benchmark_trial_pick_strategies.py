#!/usr/bin/env python3
"""Compare **PBAT converted-genome trial pickers** vs ``NM:i`` on a Bhmem/BWA BAM.

Strategies (see :func:`bhmem_equivalent_selection.pick_pbat_single_trial`):

- ``bhmem_fold`` — MAPQ/AS/NM/M ordering (usually min distance on a fixed BAM).
- ``margin_at_min_dist`` — tag-free: among global min-distance trials, maximize margin to runner-up.
- ``unique_nm_match_else_bhmem`` / ``unique_nm_match_else_nearest`` / ``nearest_nm`` — use ``NM:i``
  to disambiguate (calibration / oracle-style; not for Yap when NM is not comparable).

Example::

  python benchmark_trial_pick_strategies.py /path/to.bam /path/to/Bisulfite_Genome --max-reads 5000
"""

from __future__ import annotations

import argparse
import os
import sys

import pysam

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bhmem_equivalent_selection import (  # noqa: E402
    enumerate_pbat_single_trials,
    pick_pbat_single_trial,
)

DEFAULT_BISULFITE = (
    "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"
)

STRATEGIES = (
    "bhmem_fold",
    "margin_at_min_dist",
    "unique_nm_match_else_bhmem",
    "unique_nm_match_else_nearest",
    "nearest_nm",
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", help="BAM (primary reads with NM recommended)")
    ap.add_argument(
        "bisulfite_genome",
        nargs="?",
        default=DEFAULT_BISULFITE,
        help=f"Bismark Bisulfite_Genome dir (default: {DEFAULT_BISULFITE})",
    )
    ap.add_argument("--max-reads", type=int, default=10_000, help="Max primary mapped reads")
    args = ap.parse_args()

    ct_fa = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_fa = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for p in (ct_fa, ga_fa, args.bam):
        if not os.path.isfile(p):
            print(f"ERROR: missing file: {p}", file=sys.stderr)
            sys.exit(1)

    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)
    bam = pysam.AlignmentFile(args.bam, "rb")

    counts: dict[str, dict[str, int]] = {
        s: {"n": 0, "dist_eq_nm": 0, "label_eq_winner": 0} for s in STRATEGIES
    }
    # Reference winner: trial whose dist == NM (arbitrary if multiple — track strict label match)
    n_nm = 0
    n_eligible = 0
    n_unique_nm_trial = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("NM"):
            continue
        nm = int(read.get_tag("NM"))
        n_nm += 1
        if n_nm > args.max_reads:
            break

        trials = enumerate_pbat_single_trials(read, fa_ct, fa_ga)
        if not trials:
            continue
        n_eligible += 1
        exact_labs = [t.label for t in trials if t.dist == nm]
        if len(exact_labs) == 1:
            n_unique_nm_trial += 1
        ref_label = exact_labs[0] if len(exact_labs) == 1 else None

        for strat in STRATEGIES:
            d, lab, _reason = pick_pbat_single_trial(
                read, fa_ct, fa_ga, strategy=strat, as_per_trial=None, nm_hint=None
            )
            c = counts[strat]
            c["n"] += 1
            if d == nm:
                c["dist_eq_nm"] += 1
            if ref_label is not None and lab == ref_label:
                c["label_eq_winner"] += 1

    bam.close()
    fa_ct.close()
    fa_ga.close()

    print("BAM:", args.bam)
    print("bisulfite_genome:", args.bisulfite_genome)
    print("primary_with_NM_seen", n_nm)
    print("eligible_nonempty_trials", n_eligible)
    print("reads_with_unique_NM_matching_trial", n_unique_nm_trial)
    if n_eligible:
        print(
            "frac_eligible_unique_nm_trial",
            f"{n_unique_nm_trial / n_eligible:.6f}",
        )
    print()
    print("strategy\tn\tdist==NM\tfrac\tlabel==unique_NM_trial\tfrac_of_unique_subset")
    for strat in STRATEGIES:
        c = counts[strat]
        n = c["n"]
        if not n:
            continue
        fe = c["dist_eq_nm"] / n
        fl = (c["label_eq_winner"] / n_unique_nm_trial) if n_unique_nm_trial else float("nan")
        print(
            f"{strat}\t{n}\t{c['dist_eq_nm']}\t{fe:.6f}\t{c['label_eq_winner']}\t{fl:.6f}"
        )
    print()
    print(
        "# dist==NM: picked trial's recomputed distance matches NM:i. "
        "label==unique_NM_trial: among reads where exactly one trial has dist==NM, "
        "does the picker choose that trial's label (stronger than dist alone when min!=NM)."
    )


if __name__ == "__main__":
    main()
