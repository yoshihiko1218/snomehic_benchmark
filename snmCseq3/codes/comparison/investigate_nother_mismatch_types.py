#!/usr/bin/env python3
"""
Investigate **non–(C/T, G/A)** substitution mismatches vs **mm10 FASTA** on M/=/X columns.

``n_other`` here means: ``read_base != ref_fasta`` and the pair is **not** (C,T) and **not** (G,A).
For a sample of reads (especially read 2), prints:

  - Global counts of ``(ref_fa, read)`` pairs for read1 vs read2
  - For ``n_other`` positions: how often ``MD`` reference base **agrees with read** (true SNP in
    genome coords) vs **disagrees** (encoding / indel-adjacent / ambiguity)
  - A few **example reads** with per-position ``ref_fa``, ``read``, ``ref_md``

**Why read1 “same way” fails read2:** strand bisulfite masking only drops one conversion class per
mapping strand; mate 2 in PBAT/non-directional still shows **both** C/T and G/A–style mismatches vs
genome. That is separate from ``n_other``, which is **not** C/T or G/A at all (e.g. A↔G, T↔C,
involving **N**, etc.).

Usage:
  python investigate_nother_mismatch_types.py ref.fa alignments.bam [--mate r2] [--max-reads 5000] \\
    [--examples 5]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

def _is_ct_pair(refb: str, rb: str) -> bool:
    return refb == "C" and rb == "T"


def _is_ga_pair(refb: str, rb: str) -> bool:
    return refb == "G" and rb == "A"


def _walk_fasta_mismatches(read, fa: pysam.FastaFile) -> list[tuple[int, str, str]]:
    """List (ref_pos_0based, ref_fasta, read_base) for aligned match columns where read != ref_fasta.

    Uses ``get_aligned_pairs(matches_only=True)`` so query indices stay consistent with CIGAR/H/S.
    """
    out = []
    if read.is_unmapped or read.query_sequence is None:
        return out
    chrom = read.reference_name
    qs = read.query_sequence
    try:
        pairs = read.get_aligned_pairs(matches_only=True)
    except (ValueError, AttributeError):
        return out
    for q, rpos in pairs:
        if q is None or rpos is None:
            continue
        rb = qs[q].upper()
        refb = fa.fetch(chrom, rpos, rpos + 1).upper()
        if rb != refb:
            out.append((int(rpos), refb, rb))
    return out


def _md_ref_at_positions(read) -> dict[int, str]:
    """Map reference position -> MD-derived upper ref base (only aligned query bases)."""
    m: dict[int, str] = {}
    if not read.has_tag("MD"):
        return m
    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except (ValueError, AttributeError):
        return m
    for t in pairs:
        if len(t) != 3:
            continue
        q, rpos, refb = t
        if q is None or rpos is None or refb is None:
            continue
        m[int(rpos)] = str(refb).upper()
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description="Classify non-CT/GA substitution mismatches vs FASTA")
    ap.add_argument("reference_fasta")
    ap.add_argument("bam")
    ap.add_argument(
        "--mate",
        choices=("r1", "r2", "both"),
        default="r2",
        help="Which mate to scan (default read2)",
    )
    ap.add_argument("--max-reads", type=int, default=5000)
    ap.add_argument("--examples", type=int, default=5, help="Reads to print position-by-position")
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference_fasta)
    pair_counts: dict[str, Counter] = {
        "r1": Counter(),
        "r2": Counter(),
    }
    nother_detail: dict[str, Counter] = {
        "r1": Counter(),
        "r2": Counter(),
    }
    # For n_other: (ref_fa, read) -> count; and MD agreement
    md_read_match_on_nother = defaultdict(int)
    md_read_mismatch_on_nother = defaultdict(int)
    md_missing_on_nother = defaultdict(int)

    scored_reads: list[tuple[str, str, int, list[tuple[int, str, str, str]]]] = []
    # (qname, mate_label, n_other, list of (rpos, ref_fa, rb, ref_md or ''))

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        n = 0
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_read1:
                mate = "r1"
            elif read.is_read2:
                mate = "r2"
            else:
                continue
            if args.mate != "both" and mate != args.mate:
                continue

            mm = _walk_fasta_mismatches(read, fa)
            mdmap = _md_ref_at_positions(read)
            detail_rows: list[tuple[int, str, str, str]] = []
            n_other = 0
            for rpos, refb, rb in mm:
                key = (refb, rb)
                pair_counts[mate][key] += 1
                if _is_ct_pair(refb, rb) or _is_ga_pair(refb, rb):
                    continue
                n_other += 1
                nother_detail[mate][key] += 1
                rmd = mdmap.get(rpos, "")
                if not rmd:
                    md_missing_on_nother[mate] += 1
                elif rmd == rb:
                    md_read_match_on_nother[mate] += 1
                else:
                    md_read_mismatch_on_nother[mate] += 1
                detail_rows.append((rpos, refb, rb, rmd))

            if n_other > 0:
                scored_reads.append(
                    (read.query_name, mate, n_other, detail_rows)
                )
            n += 1
            if n >= args.max_reads:
                break

    mates_to_print = ("r1", "r2") if args.mate == "both" else (args.mate,)

    print(f"scanned_primary_reads\t{n}\tmate_filter\t{args.mate}")
    for mate in mates_to_print:
        print(f"\n## {mate}: substitution mismatch pairs vs FASTA (M/=/X only)")
        total = sum(pair_counts[mate].values())
        ct = sum(
            pair_counts[mate][k]
            for k in pair_counts[mate]
            if _is_ct_pair(k[0], k[1])
        )
        ga = sum(
            pair_counts[mate][k]
            for k in pair_counts[mate]
            if _is_ga_pair(k[0], k[1])
        )
        oth = total - ct - ga
        print(f"total_mismatch_columns\t{total}\tC/T\t{ct}\tG/A\t{ga}\tother\t{oth}")
        if oth:
            print("top (ref,read) among 'other':")
            for (a, b), c in nother_detail[mate].most_common(15):
                print(f"  {a}->{b}\t{c}")

        no = md_read_match_on_nother[mate] + md_read_mismatch_on_nother[mate] + md_missing_on_nother[mate]
        if no:
            print(
                f"\nAt 'other' columns (not ref C/read T nor ref G/read A vs mm10):"
            )
            print(
                f"  MD_ref==read\t{md_read_match_on_nother[mate]}\t"
                f"(aligner treats as match to bisulfite-aware ref; mm10 alone disagrees)"
            )
            print(
                f"  MD_ref!=read\t{md_read_mismatch_on_nother[mate]}\t"
                f"(mismatch even in MD — sequencing error / SNP / complex)"
            )
            print(
                f"  no_MD_at_pos\t{md_missing_on_nother[mate]}\t(sum={no})"
            )
            tc = nother_detail[mate].get(("T", "C"), 0)
            ag = nother_detail[mate].get(("A", "G"), 0)
            if tc or ag:
                print(
                    f"  Note: T->C columns={tc}, A->G columns={ag} — often MD has C/G vs mm10 T/A "
                    f"(conversion / opposite-strand view); not 'random' transversions."
                )

    # Examples: top n_other reads
    scored_reads.sort(key=lambda x: -x[2])
    print(f"\n## Top {args.examples} reads by n_other (vs FASTA, not C/T or G/A)")
    for read in scored_reads[: args.examples]:
        qn, mate, no, rows = read
        print(f"\n--- {qn}  {mate}  n_other={no}  (showing up to 25 positions) ---")
        for rpos, refb, rb, rmd in rows[:25]:
            print(
                f"  pos {rpos}\tref_fa={refb}\tread={rb}\tref_md={rmd or 'NA'}\t"
                f"read==mdref={rb == rmd if rmd else 'n/a'}"
            )
        if len(rows) > 25:
            print(f"  ... {len(rows) - 25} more")

    fa.close()


if __name__ == "__main__":
    main()
