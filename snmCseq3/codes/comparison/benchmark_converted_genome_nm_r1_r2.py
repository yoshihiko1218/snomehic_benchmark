#!/usr/bin/env python3
"""Benchmark converted-genome ``NM`` recompute (read1 vs read2) vs ``NM:i``.

Default bisulfite directory: ``mm10_bismark/Bisulfite_Genome`` (override with first argument).

Three estimators:

1. **md_fallback** — :func:`recompute_nm_from_converted_genomes_pbat` (trial hit ``NM`` else ``MD``).
2. **no_md_min** — :func:`recompute_nm_from_converted_genomes_pbat_no_md` (``min`` of four CT/GA×orientation trials; **no** ``MD``).
3. **no_md_pair** — :func:`recompute_nm_pair_from_converted_genomes_pbat_no_md` on fragments where both
   mates map to indexed chromosomes (joint min over Bhmem-style (CT,CT)/(GA,GA)/(CT,GA)/(GA,CT); **no** ``MD``).

Example::

  python benchmark_converted_genome_nm_r1_r2.py alignments.bam
  python benchmark_converted_genome_nm_r1_r2.py /path/to/Bisulfite_Genome alignments.bam --max-per-mate 8000
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pysam
from scipy.stats import spearmanr

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bisulfite_corrected_mismatch import (  # noqa: E402
    recompute_nm_from_converted_genomes_pbat,
    recompute_nm_from_converted_genomes_pbat_no_md,
    recompute_nm_pair_from_converted_genomes_pbat_no_md,
    recompute_nm_style_from_md,
)

DEFAULT_BISULFITE_GENOME = (
    "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"
)


def _allowed_chroms(fa_ct: pysam.FastaFile) -> set[str]:
    return {
        c.replace("_CT_converted", "")
        for c in fa_ct.references
        if c.endswith("_CT_converted")
    }


def _summarize_line(
    label: str, nm_arr: list[int], rec_arr: list[int]
) -> None:
    nm_a = np.array(nm_arr, dtype=np.int32)
    rec_a = np.array(rec_arr, dtype=np.int32)
    n = len(nm_a)
    if n == 0:
        print(f"  {label}: no reads")
        return
    eq = int(np.sum(rec_a == nm_a))
    sp, _ = spearmanr(nm_a, rec_a)
    mad = float(np.mean(np.abs(rec_a - nm_a.astype(np.float64))))
    print(
        f"  {label}: n={n}  frac(rec==NM)={eq/n:.4f}  mean|rec-NM|={mad:.4f}  Spearman={sp:.4f}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "bam_or_genome",
        nargs="?",
        default=None,
        help="If only one positional: BAM path (uses default bisulfite genome). "
        "If two positionals: Bisulfite_Genome dir then BAM.",
    )
    ap.add_argument("bam", nargs="?", default=None, help="BAM (when first arg is genome dir)")
    ap.add_argument("--max-per-mate", type=int, default=8000)
    args = ap.parse_args()

    if args.bam is not None:
        base = args.bam_or_genome.rstrip("/")
        bam_path = args.bam
    elif args.bam_or_genome is not None:
        base = DEFAULT_BISULFITE_GENOME
        bam_path = args.bam_or_genome
    else:
        ap.print_help()
        print("\nERROR: pass a BAM path.", file=sys.stderr)
        sys.exit(1)

    ct_fa = f"{base}/CT_conversion/genome_mfa.CT_conversion.fa"
    ga_fa = f"{base}/GA_conversion/genome_mfa.GA_conversion.fa"
    for p in (ct_fa, ga_fa):
        if not os.path.isfile(p):
            print(f"ERROR: missing {p}", file=sys.stderr)
            sys.exit(1)

    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)
    allowed = _allowed_chroms(fa_ct)
    print("Bisulfite genome:", base)
    print("Chromosomes in index:", sorted(allowed))
    print("BAM:", bam_path)
    print()

    # Per-mate accumulators
    r1_nm: list[int] = []
    r1_md_fb: list[int] = []
    r1_no_md: list[int] = []
    r2_nm: list[int] = []
    r2_md_fb: list[int] = []
    r2_no_md: list[int] = []

    # Pair: no_md joint (both mates on allowed chroms)
    pair_nm_sum: list[int] = []
    pair_rec_sum: list[int] = []
    pair_r1_ok = pair_r2_ok = pair_both_ok = 0

    by_name: dict[str, dict] = {}
    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_name not in allowed or not read.has_tag("NM"):
            continue

        nm = int(read.get_tag("NM"))
        rec_md = recompute_nm_from_converted_genomes_pbat(
            read, fa_ct, fa_ga, use_md_fallback=True
        )
        rec_nm = recompute_nm_from_converted_genomes_pbat_no_md(read, fa_ct, fa_ga)

        key = "r1" if read.is_read1 else "r2"
        by_name.setdefault(read.query_name, {})[key] = read

        if read.is_read1 and len(r1_nm) < args.max_per_mate:
            r1_nm.append(nm)
            r1_md_fb.append(rec_md)
            r1_no_md.append(rec_nm)
        elif read.is_read2 and len(r2_nm) < args.max_per_mate:
            r2_nm.append(nm)
            r2_md_fb.append(rec_md)
            r2_no_md.append(rec_nm)

        if len(r1_nm) >= args.max_per_mate and len(r2_nm) >= args.max_per_mate:
            break

    bam.close()

    # Second pass for pairs: scan until we have enough complete pairs on allowed chroms
    bam = pysam.AlignmentFile(bam_path, "rb")
    seen_pairs = 0
    max_pairs = 5000
    buf: dict[str, dict] = {}
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_name not in allowed or not read.has_tag("NM"):
            continue
        qn = read.query_name
        slot = buf.setdefault(qn, {})
        slot["r1" if read.is_read1 else "r2"] = read
        if "r1" in slot and "r2" in slot:
            r1, r2 = slot["r1"], slot["r2"]
            del buf[qn]
            pr = recompute_nm_pair_from_converted_genomes_pbat_no_md(r1, r2, fa_ct, fa_ga)
            if pr is None:
                continue
            d1, d2 = pr
            n1 = int(r1.get_tag("NM"))
            n2 = int(r2.get_tag("NM"))
            pair_nm_sum.append(n1 + n2)
            pair_rec_sum.append(d1 + d2)
            if d1 == n1:
                pair_r1_ok += 1
            if d2 == n2:
                pair_r2_ok += 1
            if d1 == n1 and d2 == n2:
                pair_both_ok += 1
            seen_pairs += 1
            if seen_pairs >= max_pairs:
                break
    bam.close()

    print("(1) With MD fallback when trials miss NM (tag-assisted disambiguation)")
    _summarize_line("Read 1", r1_nm, r1_md_fb)
    _summarize_line("Read 2", r2_nm, r2_md_fb)

    print()
    print("(2) No MD: min of four CT/GA × orientation trials per read")
    _summarize_line("Read 1", r1_nm, r1_no_md)
    _summarize_line("Read 2", r2_nm, r2_no_md)

    print()
    print("(3) No MD: joint pair minimization (Bhmem-style CT/GA pair types), sum(NM) vs sum(rec)")
    if pair_nm_sum:
        np_nm = np.array(pair_nm_sum, dtype=np.int32)
        np_rec = np.array(pair_rec_sum, dtype=np.int32)
        npr = len(np_nm)
        eqs = int(np.sum(np_rec == np_nm))
        sp, _ = spearmanr(np_nm, np_rec)
        mad = float(np.mean(np.abs(np_rec - np_nm.astype(np.float64))))
        print(
            f"  pairs={npr}  frac(sum_rec==sum_NM)={eqs/npr:.4f}  mean|diff|={mad:.4f}  Spearman={sp:.4f}"
        )
        print(
            f"  per-mate exact: R1 {pair_r1_ok}/{npr}  R2 {pair_r2_ok}/{npr}  both {pair_both_ok}/{npr}"
        )
    else:
        print("  no complete pairs on indexed chromosomes")

    print()
    print("Paired fragments (both mates on indexed chroms) — MD fallback per read")
    both_ok = r1_ok_r2_fail = r1_fail_r2_ok = both_fail = incomplete = 0
    for _n, d in by_name.items():
        if "r1" not in d or "r2" not in d:
            incomplete += 1
            continue
        r1, r2 = d["r1"], d["r2"]
        n1 = int(r1.get_tag("NM"))
        n2 = int(r2.get_tag("NM"))
        c1 = recompute_nm_from_converted_genomes_pbat(r1, fa_ct, fa_ga, use_md_fallback=True)
        c2 = recompute_nm_from_converted_genomes_pbat(r2, fa_ct, fa_ga, use_md_fallback=True)
        ok1, ok2 = c1 == n1, c2 == n2
        if ok1 and ok2:
            both_ok += 1
        elif ok1 and not ok2:
            r1_ok_r2_fail += 1
        elif not ok1 and ok2:
            r1_fail_r2_ok += 1
        else:
            both_fail += 1
    tot = both_ok + r1_ok_r2_fail + r1_fail_r2_ok + both_fail
    print(f"  complete pairs: {tot}  (only one mate in chr set: {incomplete})")
    if tot:
        print(f"  both rec==NM: {both_ok} ({100 * both_ok / tot:.2f}%)")
        print(f"  R1 ok, R2 miss: {r1_ok_r2_fail} ({100 * r1_ok_r2_fail / tot:.2f}%)")
        print(f"  R1 miss, R2 ok: {r1_fail_r2_ok} ({100 * r1_fail_r2_ok / tot:.2f}%)")
        print(f"  both miss: {both_fail} ({100 * both_fail / tot:.2f}%)")

    fa_ct.close()
    fa_ga.close()

    # Optional: MD vs NM (informational)
    print()
    print("Reference: frac(MD recompute == NM) on same read sample")
    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)
    bam = pysam.AlignmentFile(bam_path, "rb")
    r1_md_nm = r2_md_nm = 0
    r1_c = r2_c = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_name not in allowed or not read.has_tag("NM"):
            continue
        if not read.has_tag("MD"):
            continue
        nm = int(read.get_tag("NM"))
        mdv = recompute_nm_style_from_md(read)
        if mdv < 0:
            continue
        if read.is_read1 and r1_c < args.max_per_mate:
            r1_c += 1
            if mdv == nm:
                r1_md_nm += 1
        elif read.is_read2 and r2_c < args.max_per_mate:
            r2_c += 1
            if mdv == nm:
                r2_md_nm += 1
        if r1_c >= args.max_per_mate and r2_c >= args.max_per_mate:
            break
    bam.close()
    fa_ct.close()
    fa_ga.close()
    if r1_c:
        print(f"  Read 1: {r1_md_nm}/{r1_c} = {r1_md_nm/r1_c:.4f}")
    if r2_c:
        print(f"  Read 2: {r2_md_nm}/{r2_c} = {r2_md_nm/r2_c:.4f}")


if __name__ == "__main__":
    main()
