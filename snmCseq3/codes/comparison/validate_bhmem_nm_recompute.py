#!/usr/bin/env python3
"""
Validate NM recomputation against a Bhmem BAM.

Compares:

1. **MD + CIGAR** recompute (:func:`bisulfite_corrected_mismatch.count_nm_style_edit_distance_from_md`)
   — should match ``NM:i`` almost always on BWA/Bhmem output.
2. **Converted genome + CIGAR** (no MD): minimum over four PBAT trials, and whether ``NM`` lies in the trial set.
3. **Pair fold** (Bhmem comparator + optional enzyme bed): joint ``(d1,d2)`` vs tags.

Example::

  conda activate scnomehic
  python validate_bhmem_nm_recompute.py \\
    /path/to/SRR.bhmem.bam \\
    /path/to/Bisulfite_Genome \\
    --enzyme /path/to/dpnII.span_region.bedgraph \\
    --max-primary 25000
"""

from __future__ import annotations

import argparse
import os
import sys

import pysam

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bhmem_equivalent_selection import (  # noqa: E402
    EnzymeRegionIndex,
    recompute_nm_bhmem_style_pair_pbat_nd,
)
from bisulfite_corrected_mismatch import (  # noqa: E402
    count_nm_style_edit_distance_from_md,
    pbat_converted_genome_trial_distances,
)

DEFAULT_ENZYME = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/dpnII.span_region.bedgraph"
DEFAULT_BISULFITE = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", help="Bhmem BAM")
    ap.add_argument(
        "bisulfite_genome",
        nargs="?",
        default=DEFAULT_BISULFITE,
        help=f"Bismark Bisulfite_Genome dir (default: {DEFAULT_BISULFITE})",
    )
    ap.add_argument(
        "--enzyme",
        default=DEFAULT_ENZYME,
        help=f"Restriction BED/bedgraph for Bhmem -enzymeList (default: {DEFAULT_ENZYME})",
    )
    ap.add_argument("--no-enzyme", action="store_true", help="Skip loading enzyme index")
    ap.add_argument("--max-primary", type=int, default=30_000, help="Max primary mapped reads to scan")
    ap.add_argument("--max-pairs", type=int, default=5_000, help="Max full pairs for joint test")
    args = ap.parse_args()

    ct_fa = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_fa = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for p in (ct_fa, ga_fa, args.bam):
        if not os.path.isfile(p):
            print(f"ERROR: missing file: {p}", file=sys.stderr)
            sys.exit(1)

    enzyme: EnzymeRegionIndex | None = None
    if not args.no_enzyme and args.enzyme and os.path.isfile(args.enzyme):
        print("Loading enzyme intervals:", args.enzyme, flush=True)
        enzyme = EnzymeRegionIndex.from_bed_file(args.enzyme)
        print(f"  chromosomes_with_intervals={enzyme.n_chromosomes}", flush=True)
    elif not args.no_enzyme:
        print("WARN: enzyme file missing, continuing without enzyme:", args.enzyme, file=sys.stderr)

    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)

    n_pri = 0
    md_eq = md_miss = md_na = 0
    min_eq = min_ne = min_bad = 0
    in_set = not_in_set = 0
    pair_ok = pair_fail = pair_skip = 0

    pending: dict[str, dict] = {}

    bam = pysam.AlignmentFile(args.bam, "rb")
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("NM"):
            continue
        nm = int(read.get_tag("NM"))
        n_pri += 1
        if n_pri > args.max_primary:
            break

        if read.has_tag("MD"):
            md = count_nm_style_edit_distance_from_md(read)
            if md < 0:
                md_na += 1
            elif md == nm:
                md_eq += 1
            else:
                md_miss += 1
        else:
            md_na += 1

        trials = pbat_converted_genome_trial_distances(read, fa_ct, fa_ga)
        if trials:
            m = min(trials)
            if m == nm:
                min_eq += 1
            else:
                min_ne += 1
            if nm in trials:
                in_set += 1
            else:
                not_in_set += 1
        else:
            min_bad += 1

        if read.is_paired:
            qn = read.query_name
            slot = pending.setdefault(qn, {})
            slot["r1" if read.is_read1 else "r2"] = read
            if "r1" in slot and "r2" in slot:
                r1, r2 = slot["r1"], slot["r2"]
                del pending[qn]
                if pair_ok + pair_fail + pair_skip < args.max_pairs:
                    if (
                        r1.is_unmapped
                        or r2.is_unmapped
                        or not r1.has_tag("NM")
                        or not r2.has_tag("NM")
                    ):
                        pair_skip += 1
                    else:
                        d1, d2, _ = recompute_nm_bhmem_style_pair_pbat_nd(
                            r1, r2, fa_ct, fa_ga, enzyme=enzyme
                        )
                        n1, n2 = int(r1.get_tag("NM")), int(r2.get_tag("NM"))
                        if d1 == n1 and d2 == n2:
                            pair_ok += 1
                        else:
                            pair_fail += 1

    bam.close()
    fa_ct.close()
    fa_ga.close()

    print("BAM:", args.bam)
    print("Bisulfite genome:", args.bisulfite_genome)
    print("primary_reads_scanned", n_pri)
    print()
    print("MD_recompute vs NM:i (Bhmem should agree)")
    print(f"  equal\t{md_eq}\tnequal\t{md_miss}\tno_MD_or_fail\t{md_na}")
    if md_eq + md_miss:
        print(f"  frac_equal\t{md_eq / (md_eq + md_miss):.6f}")
    print()
    print("Converted FASTA (4 trials), no MD:")
    print(f"  min(trials)==NM\t{min_eq}\tmin!=NM\t{min_ne}\tno_trials\t{min_bad}")
    if min_eq + min_ne:
        print(f"  frac_min_eq_NM\t{min_eq / (min_eq + min_ne):.6f}")
    print(f"  NM_in_trial_set\t{in_set}\tNM_not_in_set\t{not_in_set}")
    if in_set + not_in_set:
        print(f"  frac_NM_in_set\t{in_set / (in_set + not_in_set):.6f}")
    print()
    print("Pair fold vs per-mate NM (same comparator; enzyme affects order only if loci differ):")
    print(f"  both_mates_match\t{pair_ok}\tmismatch\t{pair_fail}\tskipped\t{pair_skip}")
    print()
    print(
        "# Interpretation: MD recompute uses SEQ+MD in the BAM. Bhmem rebuilds SEQ from FASTQ "
        "via CIGAR (String2SamRecord); NM/MD still come from BWA's SAM line, so a few percent "
        "of reads can show MD_recompute != NM:i even though both are well-defined. "
        "Converted-genome 'NM_in_trial_set' checks whether NM is achievable with some CT/GA×orientation "
        "trial on the **stored** CIGAR. Pair fold uses tag AS (same for all trials); Bhmem's true "
        "pair choice uses per-candidate BWA AS, so joint (d1,d2) may not match tags without jbwa."
    )


if __name__ == "__main__":
    main()
