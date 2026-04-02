#!/usr/bin/env python3
"""
Regional **bwa mem** on bisulfite-**converted** FASTA around the BAM locus, per PBAT trial.

For each trial (CT/GA genome × converted-query orientation), builds a small reference
``[pos - pad, pos + pad)`` (clamped), ``bwa index`` + ``bwa mem`` on one read, parses
primary hit **MAPQ**, **AS:i**, **NM:i**, CIGAR-``M`` span (Bhmem-style), then folds trials
with :func:`bhmem_equivalent_selection.bhmem_prefer_second_single`.

**Local MAPQ-like:** MAPQ is computed by BWA **only against this window** (not genome-wide),
so it reflects uniqueness **within the region** — useful when you refuse to move the read
off the pipeline’s placement but still want alternative-hit signal in a neighborhood.

Compare agreement of chosen **NM** with ``NM:i`` vs path-only minimum over four converted walks
(:func:`bisulfite_corrected_mismatch.pbat_converted_genome_trial_distances`).

Requires **bwa** on PATH or ``--bwa``. Example::

  conda activate scnomehic
  python regional_bwa_trial_pick.py alignments.bam \\
    /path/to/Bisulfite_Genome \\
    --bwa /projects/.../conda/envs/scnomehic/bin/bwa \\
    --pad 8000 --max-reads 500
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pysam

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from bhmem_equivalent_selection import bhmem_prefer_second_single  # noqa: E402
from bisulfite_corrected_mismatch import (  # noqa: E402
    _pbat_converted_query_variants,
    bisulfite_converted_contig_name,
    pbat_converted_genome_trial_distances,
)

DEFAULT_BISULFITE = (
    "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"
)
DEFAULT_BWA = "/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/bwa"


def _ref_len(fa: pysam.FastaFile, contig: str) -> int:
    try:
        return int(fa.get_reference_length(contig))
    except (AttributeError, ValueError, KeyError):
        i = fa.references.index(contig)
        return int(fa.lengths[i])


def _cigar_m_len_from_string(cigar: str) -> int:
    """Bhmem / htsjdk: count per-base ``M`` only (not ``=`` / ``X``)."""
    n = 0
    for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        ln, op = int(m.group(1)), m.group(2)
        if op == "M":
            n += ln
    return n


def _parse_sam_tags(fields: list[str]) -> dict[str, str | int]:
    out: dict[str, str | int] = {}
    for f in fields:
        if f.startswith("RG:"):
            continue
        parts = f.split(":", 2)
        if len(parts) != 3:
            continue
        tag, typ, val = parts[0], parts[1], parts[2]
        if typ == "i":
            try:
                out[tag] = int(val)
            except ValueError:
                out[tag] = val
        else:
            out[tag] = val
    return out


@dataclass
class RegionalHit:
    mapq: int
    as_: int
    nm: int
    cm: int
    cigar: str


def bwa_mem_regional(
    bwa_exe: str,
    ref_seq: str,
    qname: str,
    seq: str,
    qual: str | None,
    *,
    threads: int = 1,
    expected_local_start_0: int | None = None,
    anchor_slop: int = 0,
) -> RegionalHit | None:
    """Index ``ref_seq`` as one contig ``reg``, align one read; return primary hit stats.

    If ``expected_local_start_0`` is set (BAM ``reference_start - window_start``), only accept a
    primary whose 1-based SAM POS matches within ``anchor_slop`` (keeps **same locus**).
    """
    if not ref_seq or not seq:
        return None
    if qual is None or len(qual) != len(seq):
        qual = "I" * len(seq)

    with tempfile.TemporaryDirectory(prefix="regbwa_") as td:
        ref_path = os.path.join(td, "reg.fa")
        fq_path = os.path.join(td, "read.fq")
        with open(ref_path, "w", encoding="utf-8") as f:
            f.write(">reg\n")
            f.write(ref_seq.upper())
            f.write("\n")
        r = subprocess.run(
            [bwa_exe, "index", ref_path],
            cwd=td,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if r.returncode != 0:
            return None
        with open(fq_path, "w", encoding="utf-8") as f:
            f.write(f"@{qname}\n{seq.upper()}\n+\n{qual}\n")
        r2 = subprocess.run(
            [
                bwa_exe,
                "mem",
                "-t",
                str(threads),
                ref_path,
                fq_path,
            ],
            cwd=td,
            capture_output=True,
            text=True,
            timeout=180,
        )
        if r2.returncode != 0:
            return None

        best: RegionalHit | None = None
        best_d = 10**9
        for line in r2.stdout.splitlines():
            if line.startswith("@") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 11:
                continue
            fl = int(p[1])
            if fl & 4:
                continue
            if fl & 0x100:
                continue
            if fl & 0x800:
                continue
            pos1 = int(p[3])
            pos0 = pos1 - 1
            if expected_local_start_0 is not None:
                d = abs(pos0 - expected_local_start_0)
                if d > anchor_slop:
                    continue
            else:
                d = 0
            mq = int(p[4])
            cigar = p[5]
            tags = _parse_sam_tags(p[11:])
            if "AS" not in tags or "NM" not in tags:
                continue
            as_ = int(tags["AS"])
            nm = int(tags["NM"])
            cm = _cigar_m_len_from_string(cigar)
            hit = RegionalHit(mapq=mq, as_=as_, nm=nm, cm=cm, cigar=cigar)
            if expected_local_start_0 is None:
                return hit
            if best is None or d < best_d or (
                d == best_d
                and bhmem_prefer_second_single(
                    best.mapq, best.as_, best.nm, best.cm, mq, as_, nm, cm
                )
            ):
                best, best_d = hit, d
        return best


def fold_trials_bhmem(
    trials: list[tuple[str, str, RegionalHit]],
) -> tuple[RegionalHit, str, str] | None:
    """Return (winning hit, genome, qtag)."""
    if not trials:
        return None
    best_g, best_t, best_h = trials[0]
    bmq, bas, bnm, bcm = best_h.mapq, best_h.as_, best_h.nm, best_h.cm
    for g, t, h in trials[1:]:
        mq, as_, nm, cm = h.mapq, h.as_, h.nm, h.cm
        if bhmem_prefer_second_single(bmq, bas, bnm, bcm, mq, as_, nm, cm):
            best_g, best_t, best_h = g, t, h
            bmq, bas, bnm, bcm = mq, as_, nm, cm
    return (best_h, best_g, best_t)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("bisulfite_genome", nargs="?", default=DEFAULT_BISULFITE)
    ap.add_argument("--bwa", default="", help=f"default tries {DEFAULT_BWA}")
    ap.add_argument("--pad", type=int, default=8000, help="bp padding each side of alignment span")
    ap.add_argument(
        "--anchor-slop",
        type=int,
        default=-1,
        help="If >= 0, keep **same locus**: only accept hits whose start is within this many bp "
        "of the BAM start in window coordinates. -1 = use BWA primary as-is (can move within window).",
    )
    ap.add_argument("--max-reads", type=int, default=400)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    bwa = args.bwa or DEFAULT_BWA
    if not shutil.which(bwa) and not os.path.isfile(bwa):
        print("ERROR: bwa not found; pass --bwa path", file=sys.stderr)
        sys.exit(1)

    ct_fa = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_fa = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for p in (ct_fa, ga_fa, args.bam):
        if not os.path.isfile(p):
            print(f"ERROR: missing {p}", file=sys.stderr)
            sys.exit(1)

    fa_ct = pysam.FastaFile(ct_fa)
    fa_ga = pysam.FastaFile(ga_fa)
    bam = pysam.AlignmentFile(args.bam, "rb")

    n = reg_eq = reg_ne = pathmin_eq = pathmin_ne = path_in_set = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None or not read.has_tag("NM"):
            continue

        nm_tag = int(read.get_tag("NM"))
        qs = read.query_sequence
        if read.query_qualities is not None and len(read.query_qualities) == len(qs):
            qual = "".join(chr(q + 33) for q in read.query_qualities)
        else:
            qual = None

        rs = int(read.reference_start)
        re_ = int(read.reference_end)
        trials_hits: list[tuple[str, str, RegionalHit]] = []

        for conv, fa in ("CT", fa_ct), ("GA", fa_ga):
            cname = bisulfite_converted_contig_name(fa, read.reference_name, conv)
            if cname is None:
                continue
            rlen = _ref_len(fa, cname)
            w0 = max(0, rs - args.pad)
            w1 = min(rlen, re_ + args.pad)
            if w1 <= w0:
                continue
            ref_slice = fa.fetch(cname, w0, w1)
            exp_local = (rs - w0) if args.anchor_slop >= 0 else None
            slop = int(args.anchor_slop) if args.anchor_slop >= 0 else 0
            for qtag, qconv in _pbat_converted_query_variants(read):
                if len(qconv) != len(qs):
                    continue
                hit = bwa_mem_regional(
                    bwa,
                    ref_slice,
                    read.query_name.replace("\t", "_")[:200],
                    qconv,
                    qual,
                    threads=args.threads,
                    expected_local_start_0=exp_local,
                    anchor_slop=slop,
                )
                if hit is None:
                    continue
                trials_hits.append((conv, qtag, hit))

        trials_dists = pbat_converted_genome_trial_distances(read, fa_ct, fa_ga)
        path_min = min(trials_dists) if trials_dists else -1
        in_set = nm_tag in trials_dists if trials_dists else False

        picked = fold_trials_bhmem(trials_hits)
        if picked is None:
            continue
        win_h, _g, _t = picked

        n += 1
        if win_h.nm == nm_tag:
            reg_eq += 1
        else:
            reg_ne += 1
        if path_min == nm_tag:
            pathmin_eq += 1
        else:
            pathmin_ne += 1
        if in_set:
            path_in_set += 1

        if n >= args.max_reads:
            break

    bam.close()
    fa_ct.close()
    fa_ga.close()

    print("BAM:", args.bam)
    print(
        "pad_bp",
        args.pad,
        "anchor_slop",
        args.anchor_slop,
        "max_reads",
        n,
    )
    print("regional_bwa_bhmem_pick  NM==tag", reg_eq, " NM!=tag", reg_ne)
    if n:
        print("  frac_reg_eq", reg_eq / n)
    print("path_min_dist        NM==tag", pathmin_eq, " NM!=tag", pathmin_ne)
    if n:
        print("  frac_pathmin_eq", pathmin_eq / n)
    print("path_NM_in_4trial_set", path_in_set, " of", n)
    if n:
        print("  frac_in_set", path_in_set / n)


if __name__ == "__main__":
    main()
