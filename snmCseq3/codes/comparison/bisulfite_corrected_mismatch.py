#!/usr/bin/env python3
"""
Compute bisulfite-aware "corrected" mismatch counts from aligned BAM + reference FASTA.

**Matching ``NM:i`` on bhmem, then the same logic on yap**

If the goal is to **recompute a value as close as possible to each BAM's ``NM:i`` tag** using
**one identical implementation** for both pipelines, use :func:`recompute_nm_style_from_md`
(:func:`count_nm_style_edit_distance_from_md`). It decodes substitution columns from ``MD:Z`` and
adds CIGAR insertion/deletion lengths — no mm10 FASTA walk, no yap-vs-bhmem branches. On bhmem this
typically tracks ``NM`` very closely; on Bismark/yap it usually matches ``NM`` exactly because
``NM`` and ``MD`` are generated together. Cross-pipeline **comparison** is then between two
aligner-native edit distances (not the same as a single raw-genomic definition).

You can also recompute in **BWA index space** using **bisulfite-converted** genome FASTAs
(``*_CT_converted`` / ``*_GA_converted``) plus a **converted** query string; see
:func:`count_nm_style_edit_distance_converted_explicit` and
:func:`recompute_nm_from_converted_genomes_pbat`. That matches ``bwa_gen_cigar2`` only when the
CT vs GA genome and query-orientation match the winning alignment; ``NM:i`` (or ``MD``) disambiguates.

To pick the winning trial like **Bhmem** (MAPQ / ``AS`` / ``NM`` / CIGAR-``M`` tie-break) instead of
``min(trials)``, use ``bhmem_equivalent_selection`` (``recompute_nm_bhmem_style_single_pbat``,
``recompute_nm_bhmem_style_pair_pbat_nd``).

Genomic FASTA walks with bisulfite masking (``count_nm_style_edit_distance(...,
bisulfite_correct=True, ...)``) remain **much weaker** proxies for bhmem ``NM`` on mate 2 / pooled
libraries than ``MD``-based recompute; see ``report_nm_recompute_by_mate.py``.

Why NM differs (yap/Bismark vs bhmem):
  - Bismark aligns reads to bisulfite-converted Bowtie indices but writes SAM against the
    *genomic* reference. NM/MD count differences vs the *unconverted* genome, so
    unmethylated C -> T in the read appears as mismatches (often most of NM).
  - Bhmem/BWA-style bisulfite aligners index converted genomes and report alignments where
    those conversions are already "matches" to the index, so NM stays low.

Corrected mismatch heuristic (strand-aware, assumes full bisulfite conversion at unmethylated C):
  - Read mapped to forward genomic strand (FLAG & 16 == 0): do not count ref C vs read T.
  - Read mapped to reverse strand (FLAG & 16): do not count ref G vs read A (complement of C->T).

For PBAT / non-directional paired-end (e.g. bhmem ``-pbat -nonDirectional``), mate 2 often disagrees
with ``NM`` under strand-only masking. ``bisulfite_read2_mode="pbat_read2"`` masks **both** (C,T) and
(G,A) on read 2 (strand-agnostic). ``bisulfite_read2_mode="pbat_r2_fwd_ga_rev_ct"`` (default for
:func:`count_non_md_cross_pipeline_pbat_nd_edit_distance`) uses read-2 **mapping strand**: skip G/A
and A/G on forward read 2, skip C/T and T/C on reverse read 2 — better Spearman vs bhmem ``NM`` on
mate 2 in benchmarks while keeping read 1 on the strand rule.

This ignores methylation (true mC still looks like C in ref vs T in read and would be
subtracted too — for strict methylation-aware correction use methylation callers / XM tags).

Usage:
  python bisulfite_corrected_mismatch.py ref.fa alignments.bam [--max-reads N]
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass

try:
    import pysam
except ImportError:
    print("ERROR: pip install pysam", file=sys.stderr)
    sys.exit(1)

# CIGAR: only these consume both query and reference (substitution / match columns).
# Excludes soft clip (4), hard clip (5), insertion (1), deletion (2), skip (3), etc.
_CIGAR_MATCH_MISMATCH_OPS = frozenset(
    (
        0,
        7,
        8,
    )  # M, =, X (pysam / BAM)
)

_DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _reverse_complement_dna(seq: str) -> str:
    return seq.upper().translate(_DNA_COMPLEMENT)[::-1]


def bisulfite_converted_contig_name(
    fasta: pysam.FastaFile, chrom: str, conversion: str
) -> str | None:
    """Resolve Bismark-style bisulfite FASTA contig (e.g. ``chr1`` → ``chr1_CT_converted``).

    ``conversion`` must be ``\"CT\"`` or ``\"GA\"``. Returns ``None`` if no matching contig exists.
    """
    c = conversion.upper()
    if c == "CT":
        suf = "_CT_converted"
    elif c == "GA":
        suf = "_GA_converted"
    else:
        raise ValueError('conversion must be "CT" or "GA"')
    cand = chrom + suf
    if cand in fasta.references:
        return cand
    if chrom in fasta.references:
        return chrom
    return None


def _pbat_converted_query_variants(read) -> list[tuple[str, str]]:
    """Two converted-query strings in **SAM ``query_sequence`` order** (Bhmem ``-pbat``, not ``-snm3c``).

    BWA compared mate 1 after **G→A** and mate 2 after **C→T** on the sequencing read. Depending on
    strand and Bhmem's SAM encoding, either a direct replace or ``RC(convert(RC(SEQ)))`` matches the
    bytes BWA used; both are tried by :func:`recompute_nm_from_converted_genomes_pbat`.
    """
    qs = read.query_sequence.upper()
    if read.is_read2:

        def conv(s: str) -> str:
            return s.replace("C", "T")

    else:

        def conv(s: str) -> str:
            return s.replace("G", "A")

    v1 = conv(qs)
    v2 = _reverse_complement_dna(conv(_reverse_complement_dna(qs)))
    return [("seq_conv", v1), ("rc_conv_rc", v2)]


def count_mismatches_corrected_aligned_pairs(read, fasta):
    """Substitution columns via pysam aligned pairs (matches_only=True).

    Excludes pure indels and bases outside the aligned path (e.g. soft clip).
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0, 0, 0

    chrom = read.reference_name
    raw_mm = 0
    corr_mm = 0
    bs_ign = 0
    rev = read.is_reverse

    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if qpos is None or rpos is None:
            continue
        rb = read.query_sequence[qpos].upper()
        refb = fasta.fetch(chrom, rpos, rpos + 1).upper()
        if rb == refb:
            continue
        raw_mm += 1
        if not rev:
            if refb == "C" and rb == "T":
                bs_ign += 1
                continue
        else:
            if refb == "G" and rb == "A":
                bs_ign += 1
                continue
        corr_mm += 1

    return raw_mm, corr_mm, bs_ign


def count_mismatches_corrected_cigar_mx(read, fasta):
    """Bisulfite-corrected mismatch using **only** CIGAR M / = / X blocks.

    Walks the CIGAR string explicitly: soft-clipped (S/H) and indels (I/D/N) do not
    contribute substitution counts. This matches the usual definition of “aligned
    match/mismatch columns” and is stricter than counting NM (which includes indels).

    Returns (raw_mismatch_count, corrected_mismatch_count, bisulfite_ignored).
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0, 0, 0
    ct = read.cigartuples
    if not ct:
        return 0, 0, 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return count_mismatches_corrected_aligned_pairs(read, fasta)

    chrom = read.reference_name
    qs = read.query_sequence
    rev = read.is_reverse
    raw_mm = 0
    corr_mm = 0
    bs_ign = 0

    q = 0
    r = read.reference_start

    for op, length in ct:
        if op in _CIGAR_MATCH_MISMATCH_OPS:
            for _ in range(length):
                rb = qs[q].upper()
                refb = fasta.fetch(chrom, r, r + 1).upper()
                if rb != refb:
                    raw_mm += 1
                    if not rev:
                        if refb == "C" and rb == "T":
                            bs_ign += 1
                        else:
                            corr_mm += 1
                    else:
                        if refb == "G" and rb == "A":
                            bs_ign += 1
                        else:
                            corr_mm += 1
                q += 1
                r += 1
        elif op == 1:
            q += length
        elif op in (2, 3):
            r += length
        elif op == 4:
            q += length
        elif op == 5:
            pass

    return raw_mm, corr_mm, bs_ign


def count_mismatches_corrected(read, fasta):
    """Default: aligned-pairs method (same as original script)."""
    return count_mismatches_corrected_aligned_pairs(read, fasta)


def _bisulfite_skip_substitution(
    refb: str,
    rb: str,
    rev: bool,
    *,
    bisulfite_read2_mode: str,
    is_read2: bool,
) -> bool:
    """Return True if this substitution should not count toward edit distance (bisulfite-style).

    **Read 1 / unpaired:** strand-only masking (forward → skip C/T; reverse → skip G/A).

    **Read 2:** depends on ``bisulfite_read2_mode`` when ``is_read2``:

    - ``"strand"``: same strand rule as read 1 (often poor vs bhmem ``NM`` on mate 2).
    - ``"pbat_read2"``: skip both (C,T) and (G,A) vs genome (strand-agnostic).
    - ``"pbat_r2_fwd_ga_rev_ct"``: forward read 2 skip G/A and A/G; reverse read 2 skip C/T and T/C
      (heuristic tuned vs bhmem ``NM`` without ``MD``; see ``benchmark_bisulfite_nm_hypotheses.py``).

    Mate-2 ``NM`` still only approximates a raw FASTA walk; use ``MD``-based recompute for exact tag
    agreement.
    """
    if bisulfite_read2_mode == "pbat_read2" and is_read2:
        return (refb == "C" and rb == "T") or (refb == "G" and rb == "A")
    if bisulfite_read2_mode == "pbat_r2_fwd_ga_rev_ct" and is_read2:
        if not rev:
            return (refb == "G" and rb == "A") or (refb == "A" and rb == "G")
        return (refb == "C" and rb == "T") or (refb == "T" and rb == "C")
    if not rev:
        return refb == "C" and rb == "T"
    return refb == "G" and rb == "A"


def _bisulfite_skip_cross_pipeline_unified(refb: str, rb: str, read) -> bool:
    """Conversion skip rule for **one** non-MD metric shared across yap (Bismark) and bhmem.

    - If ``XR:Z`` is present and is ``CT`` or ``GA`` (Bismark), use that conversion class only.
    - Otherwise treat as bhmem-style: strand on read 1; on read 2 use
      ``bisulfite_read2_mode="pbat_r2_fwd_ga_rev_ct"`` (no ``XR`` on bhmem).

    This does **not** reproduce ``NM:i`` on yap BAMs (Bismark ``NM`` is uncorrected vs genome).
    On bhmem, ``NM`` is bisulfite-aware; this walk only **approximates** it (read1 is typically
    closest). Use ``count_nm_style_edit_distance_from_md`` when the goal is to match ``NM`` exactly.
    """
    if read.has_tag("XR"):
        xr = str(read.get_tag("XR")).upper()
        if xr == "CT":
            return refb == "C" and rb == "T"
        if xr == "GA":
            return refb == "G" and rb == "A"
    rev = read.is_reverse
    is_r2 = bool(read.is_paired and read.is_read2)
    mode = "pbat_r2_fwd_ga_rev_ct" if is_r2 else "strand"
    return _bisulfite_skip_substitution(
        refb, rb, rev, bisulfite_read2_mode=mode, is_read2=is_r2
    )


def count_cross_pipeline_comparable_edit_distance(read, fasta) -> int:
    """Non-MD edit distance: **genomic** FASTA + CIGAR, with unified bisulfite masking (see above).

    Returns substitution mismatches (after masking) + insertion + deletion/skip lengths, same as
    ``count_nm_style_edit_distance`` indel handling. Soft/hard clip excluded.

    Intended so **yap** and **bhmem** can be compared on the **same definition** without using
    ``MD:Z``. Expect **different** concordance with each pipeline's ``NM:i`` tag.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0
    ct = read.cigartuples
    if not ct:
        return 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return -1

    chrom = read.reference_name
    qs = read.query_sequence
    subs = 0
    indel_bases = 0
    q = 0
    r = read.reference_start

    for op, length in ct:
        if op in _CIGAR_MATCH_MISMATCH_OPS:
            for _ in range(length):
                rb = qs[q].upper()
                refb = fasta.fetch(chrom, r, r + 1).upper()
                if rb != refb:
                    if not _bisulfite_skip_cross_pipeline_unified(refb, rb, read):
                        subs += 1
                q += 1
                r += 1
        elif op == 1:
            indel_bases += length
            q += length
        elif op in (2, 3):
            indel_bases += length
            r += length
        elif op == 4:
            q += length
        elif op == 5:
            pass

    return subs + indel_bases


def count_nm_style_edit_distance(
    read,
    fasta,
    bisulfite_correct: bool = False,
    *,
    bisulfite_read2_mode: str = "strand",
) -> int:
    """Recompute SAM-style edit distance: substitution differences + insertion + deletion bases.

    Walks CIGAR: M/=/X columns compare read vs reference; I adds ``length``, D/N add ``length``.
    Soft/Hard clip excluded from the distance (same idea as ``NM``).

    If ``bisulfite_correct`` is True, substitution positions that look like unmethylated
    conversion do **not** add +1 to the substitution part (indel lengths unchanged).

    Default (``bisulfite_read2_mode="strand"``): ref C vs read T when not ``is_reverse``;
    ref G vs read A when ``is_reverse`` (complement of C→T on the other strand).

    ``bisulfite_read2_mode="pbat_read2"``: for **paired** reads with ``FLAG`` read2 (mate 2), treat
    **both** (C,T) and (G,A) substitution mismatches as conversion (strand-agnostic). Read 1 and
    unpaired reads still use the strand rule.

    ``bisulfite_read2_mode="pbat_r2_fwd_ga_rev_ct"``: read 2 only — forward mate skips G/A and A/G;
    reverse mate skips C/T and T/C; read 1 / unpaired use the strand rule. Tends to track bhmem
    ``NM`` on read 2 better than ``pbat_read2`` in Spearman benchmarks (still heuristic).

    Returns a single integer comparable to ``NM:i`` only when ``bisulfite_correct`` is
    False and the walk matches the aligner's reference bases.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0
    ct = read.cigartuples
    if not ct:
        return 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return -1

    if bisulfite_read2_mode not in (
        "strand",
        "pbat_read2",
        "pbat_r2_fwd_ga_rev_ct",
    ):
        raise ValueError(
            "bisulfite_read2_mode must be 'strand', 'pbat_read2', or "
            "'pbat_r2_fwd_ga_rev_ct'"
        )

    chrom = read.reference_name
    qs = read.query_sequence
    rev = read.is_reverse
    is_read2 = bool(read.is_paired and read.is_read2)
    subs = 0
    indel_bases = 0

    q = 0
    r = read.reference_start

    for op, length in ct:
        if op in _CIGAR_MATCH_MISMATCH_OPS:
            for _ in range(length):
                rb = qs[q].upper()
                refb = fasta.fetch(chrom, r, r + 1).upper()
                if rb != refb:
                    if bisulfite_correct:
                        if not _bisulfite_skip_substitution(
                            refb,
                            rb,
                            rev,
                            bisulfite_read2_mode=bisulfite_read2_mode,
                            is_read2=is_read2,
                        ):
                            subs += 1
                    else:
                        subs += 1
                q += 1
                r += 1
        elif op == 1:
            indel_bases += length
            q += length
        elif op in (2, 3):
            indel_bases += length
            r += length
        elif op == 4:
            q += length
        elif op == 5:
            pass

    return subs + indel_bases


def count_non_md_cross_pipeline_pbat_nd_edit_distance(read, fasta) -> int:
    """**Non-MD** edit distance (mm10 + CIGAR), **same call** on bhmem and yap for fair comparison.

    Tuned for libraries aligned like Bhmem with ``-pbat`` and ``-nonDirectional`` (see
    ``snmCseq3/codes/02.alignment.sh``):

    - **Read 1 / unpaired:** strand bisulfite mask (ignore ref C/read T forward; ref G/read A reverse).
    - **Read 2:** ``pbat_r2_fwd_ga_rev_ct`` — forward mate skips G/A and A/G; reverse mate skips C/T
      and T/C (better Spearman vs bhmem ``NM`` on mate 2 than symmetric ``pbat_read2`` in benchmarks).

    No ``MD:Z`` required. On yap/Bismark, if ``XR:Z`` is present, :func:`count_cross_pipeline_comparable_edit_distance`
    uses CT/GA from that tag instead; this function stays **XR-free** for a single CIGAR+FASTA rule.

    **Bhmem:** ``NM`` agreement is **approximate** on mate 2 (FASTA vs index); read 1 remains the closer mate.

    **Yap/Bismark:** will **not** match recorded ``NM:i`` (uncorrected vs genome). Use for a shared
    definition when comparing pipelines, not for recovering yap tags.

    Returns ``-1`` if CIGAR contains unsupported ops (same as ``count_nm_style_edit_distance``).
    """
    return count_nm_style_edit_distance(
        read,
        fasta,
        True,
        bisulfite_read2_mode="pbat_r2_fwd_ga_rev_ct",
    )


def count_nm_style_edit_distance_converted_explicit(
    read,
    conv_fasta: pysam.FastaFile,
    *,
    ref_contig: str,
    query_converted: str,
) -> int:
    """SAM-style edit distance: **converted read** vs **converted reference** + CIGAR indels.

    This is the same object BWA scores in ``bwa_gen_cigar2`` for a bisulfite index hit:
    substitution columns compare ``query_converted[qpos]`` to ``conv_fasta`` at ``rpos`` (genomic
    coordinates, same numbering as the unconverted assembly), plus insertion and deletion/skip
    lengths from the CIGAR.

    ``query_converted`` must be the same length as ``read.query_sequence`` and in **SAM query
    order**. Build it with mate-specific conversion (PBAT: R1 G→A, R2 C→T) and the correct
    orientation; :func:`recompute_nm_from_converted_genomes_pbat` tries two orientations.

    Returns ``-1`` if CIGAR is unsupported or contig/length checks fail.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0
    qs = read.query_sequence
    if len(query_converted) != len(qs):
        return -1
    qc = query_converted.upper()
    ct = read.cigartuples
    if not ct:
        return 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return -1
    if ref_contig not in conv_fasta.references:
        return -1

    subs = 0
    try:
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if qpos is None or rpos is None:
                continue
            rb = qc[qpos]
            refb = conv_fasta.fetch(ref_contig, rpos, rpos + 1).upper()
            if rb != refb:
                subs += 1
    except (ValueError, KeyError):
        return -1

    indel_bases = 0
    for op, length in ct:
        if op == 1:
            indel_bases += length
        elif op in (2, 3):
            indel_bases += length

    return subs + indel_bases


def pbat_converted_genome_trial_distances(
    read,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> list[int]:
    """Four distances: CT/GA genome × two PBAT converted-query shapes (Bhmem ``-pbat``, not ``-snm3c``).

    Empty list if contigs are missing or all trials fail. Does **not** use ``MD``.
    """
    if read.is_unmapped or read.query_sequence is None:
        return []

    c_ct = bisulfite_converted_contig_name(fasta_ct, read.reference_name, "CT")
    c_ga = bisulfite_converted_contig_name(fasta_ga, read.reference_name, "GA")
    if c_ct is None or c_ga is None:
        return []

    trials: list[int] = []
    for fa, cname in (fasta_ct, c_ct), (fasta_ga, c_ga):
        for _tag, qconv in _pbat_converted_query_variants(read):
            d = count_nm_style_edit_distance_converted_explicit(
                read, fa, ref_contig=cname, query_converted=qconv
            )
            if d >= 0:
                trials.append(d)
    return trials


def recompute_nm_from_converted_genomes_pbat_no_md(
    read,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> int:
    """Recover ``NM``-style distance **without** reading ``MD:Z``.

    **What Bhmem actually does (alignment time):** for each mate it runs **two** BWA index passes
    (CT vs GA, with PBAT naming swap as in ``Bhmem.java``), and within **each** pass it keeps the
    best **SAM line** using ``comparingSamRecord``: higher MAPQ, then higher ``AS``, then **lower**
    ``NM``, then more CIGAR ``M`` span. It does **not** score candidates by recomputing edit
    distance on a single fixed CIGAR across four converted-query trials.

    **What this function does (post hoc on an existing BAM):** it builds the same **four**
    converted-genome trials as :func:`pbat_converted_genome_trial_distances`, then applies
    :func:`bhmem_equivalent_selection.recompute_nm_bhmem_style_single_pbat` — i.e. Bhmem’s
    **single-record** ordering on those **synthetic** trials. On a **fixed** alignment, every
    trial shares the same MAPQ, CIGAR, and usually the same ``AS:i`` tag, so that ordering
    **collapses to choosing the minimum recomputed distance** (equivalently the old
    ``min(trials)``). The returned integer matches ``min(trials)``; the fold only adds Bhmem’s
    tie-breakers if those fields ever differ per trial (e.g. future jbwa ``as_per_trial``).

    **Not identical to full Bhmem** when BWA would have emitted **different** CIGARs/positions per
    index pass; reproducing that requires realigning (e.g. jbwa) and the same pair logic for PE.

    Returns ``-1`` if no valid trial. Unmapped reads return ``0``.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0
    from bhmem_equivalent_selection import recompute_nm_bhmem_style_single_pbat

    d, _ = recompute_nm_bhmem_style_single_pbat(read, fasta_ct, fasta_ga, as_per_trial=None)
    return d


def recompute_nm_pair_from_converted_genomes_pbat_no_md(
    read1,
    read2,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> tuple[int, int] | None:
    """Joint **no-MD** heuristic: minimize ``d1 + d2`` over conversion assignments.

    Enumerates **(CT,CT), (GA,GA), (CT,GA), (GA,CT)** × PBAT query-orientation variants per mate.
    Chooses the assignment with **smallest** ``d1 + d2``. This is **not** Bhmem’s
    ``comparingSamRecordPbat`` (sum MAPQ, enzyme, sum ``AS``, sum ``NM``, …). For that, use
    :func:`bhmem_equivalent_selection.recompute_nm_bhmem_style_pair_pbat_nd`.

    Returns ``None`` if mates are on chromosomes missing from the bisulfite FASTAs or no valid
    assignment exists.
    """
    if (
        read1.is_unmapped
        or read2.is_unmapped
        or read1.query_sequence is None
        or read2.query_sequence is None
    ):
        return None

    pair_index = (("CT", "CT"), ("GA", "GA"), ("CT", "GA"), ("GA", "CT"))
    best: tuple[int, int, int] | None = None  # (sum, d1, d2)

    for idx1, idx2 in pair_index:
        fa1 = fasta_ct if idx1 == "CT" else fasta_ga
        fa2 = fasta_ct if idx2 == "CT" else fasta_ga
        c1 = bisulfite_converted_contig_name(fa1, read1.reference_name, idx1)
        c2 = bisulfite_converted_contig_name(fa2, read2.reference_name, idx2)
        if c1 is None or c2 is None:
            continue
        for _t1, q1 in _pbat_converted_query_variants(read1):
            d1 = count_nm_style_edit_distance_converted_explicit(
                read1, fa1, ref_contig=c1, query_converted=q1
            )
            if d1 < 0:
                continue
            for _t2, q2 in _pbat_converted_query_variants(read2):
                d2 = count_nm_style_edit_distance_converted_explicit(
                    read2, fa2, ref_contig=c2, query_converted=q2
                )
                if d2 < 0:
                    continue
                s = d1 + d2
                if best is None or s < best[0]:
                    best = (s, d1, d2)

    if best is None:
        return None
    return (best[1], best[2])


def recompute_nm_from_converted_genomes_pbat(
    read,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
    *,
    use_md_fallback: bool = True,
) -> int:
    """Recompute ``NM``-style distance using **CT + GA converted** FASTAs (BWA bisulfite index space).

    For Bhmem **-pbat** (ignore ``-snm3c``): tries **CT** and **GA** converted references and two
    converted-query orientations per :func:`_pbat_converted_query_variants`. If ``NM:i`` is present
    and **any** trial equals it, returns that value (they should agree with the tag). If
    ``use_md_fallback`` is True (default) and ``MD:Z`` exists, falls back to
    :func:`recompute_nm_style_from_md` when no trial equals ``NM``. If ``use_md_fallback`` is False,
    uses :func:`recompute_nm_from_converted_genomes_pbat_no_md` (``min`` of trials) when no trial
    equals ``NM``. Otherwise returns ``-1``.

    Requires bisulfite genome contigs such as ``chr1_CT_converted`` / ``chr1_GA_converted`` (see
    :func:`bisulfite_converted_contig_name`) covering the BAM's ``reference_name``.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0

    trials = pbat_converted_genome_trial_distances(read, fasta_ct, fasta_ga)
    if not trials:
        if use_md_fallback and read.has_tag("MD"):
            return recompute_nm_style_from_md(read)
        return -1

    if read.has_tag("NM"):
        nm = int(read.get_tag("NM"))
        if nm in trials:
            return nm

    if use_md_fallback and read.has_tag("MD"):
        return recompute_nm_style_from_md(read)

    return recompute_nm_from_converted_genomes_pbat_no_md(read, fasta_ct, fasta_ga)


def count_nm_style_edit_distance_from_md(read) -> int:
    """SAM-style edit distance using **reference bases from the MD tag** (via ``get_aligned_pairs(with_seq=True)``).

    Counts substitution mismatches where read base differs from the **MD-encoded** reference base at
    each aligned column, plus insertion and deletion lengths (same CIGAR indel handling as
    ``count_nm_style_edit_distance``).

    For BWA/bhmem alignments with ``MD:Z``, this typically **matches ``NM:i``** (often >99.9% exact)
    because ``NM`` and ``MD`` are defined from the **same** effective reference the aligner used
    (bisulfite-aware: e.g. C→T conversion columns match the index, not raw genomic ``C`` vs read ``T``
    from ``.fa``).

    A **genomic FASTA** walk (``count_nm_style_edit_distance(..., bisulfite_correct=False)``) can
    **disagree** strongly with ``NM`` at those sites; strand bisulfite masking helps but does not
    reproduce ``NM`` as closely as this MD-based recompute.

    Returns ``-1`` if ``MD`` is missing or pairs cannot be built.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0
    if not read.has_tag("MD"):
        return -1
    ct = read.cigartuples
    if not ct:
        return 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return -1

    qs = read.query_sequence
    subs = 0
    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except (ValueError, AttributeError):
        return -1

    for t in pairs:
        if len(t) != 3:
            continue
        q, rpos, refb = t
        if q is None or rpos is None or refb is None:
            continue
        if qs[q].upper() != refb.upper():
            subs += 1

    indel_bases = 0
    for op, length in ct:
        if op == 1:
            indel_bases += length
        elif op in (2, 3):
            indel_bases += length

    return subs + indel_bases


def recompute_nm_style_from_md(read) -> int:
    """Recompute SAM edit distance from ``MD:Z`` + CIGAR indels (alias of ``count_nm_style_edit_distance_from_md``).

    **Use this same function on bhmem and on yap** when you want:

    1. A recomputation that **agrees closely with bhmem ``NM:i``** (because ``MD`` encodes the same
       effective reference the bisulfite-aware aligner used).
    2. The **identical code path** on yap/Bismark BAMs (typically **exact** agreement with yap
       ``NM:i``).

    Returns ``-1`` if ``MD`` is missing or aligned pairs cannot be built.
    """
    return count_nm_style_edit_distance_from_md(read)


def nm_style_distance_breakdown(read, fasta):
    """Decompose genomic edit distance along the alignment (same walk as ``count_nm_style_edit_distance``).

    Returns ``(n_sub_ct, n_sub_other, indel_bases)`` where:

    - ``n_sub_ct``: substitution columns that match **unmethylated conversion** vs genome
      (ref C vs read T on forward; ref G vs read A on reverse).
    - ``n_sub_other``: all other substitution mismatches.
    - ``indel_bases``: sum of I lengths + D/N lengths.

    ``n_sub_ct + n_sub_other + indel_bases`` equals ``count_nm_style_edit_distance(..., False)``.
    """
    if read.is_unmapped or read.query_sequence is None:
        return 0, 0, 0
    ct = read.cigartuples
    if not ct:
        return 0, 0, 0
    for op, _ in ct:
        if op not in (0, 1, 2, 3, 4, 5, 7, 8):
            return -1, -1, -1

    chrom = read.reference_name
    qs = read.query_sequence
    rev = read.is_reverse
    n_ct = n_other = 0
    indel_bases = 0
    q = 0
    r = read.reference_start

    for op, length in ct:
        if op in _CIGAR_MATCH_MISMATCH_OPS:
            for _ in range(length):
                rb = qs[q].upper()
                refb = fasta.fetch(chrom, r, r + 1).upper()
                if rb != refb:
                    is_ct = (not rev and refb == "C" and rb == "T") or (
                        rev and refb == "G" and rb == "A"
                    )
                    if is_ct:
                        n_ct += 1
                    else:
                        n_other += 1
                q += 1
                r += 1
        elif op == 1:
            indel_bases += length
            q += length
        elif op in (2, 3):
            indel_bases += length
            r += length
        elif op == 4:
            q += length
        elif op == 5:
            pass

    return n_ct, n_other, indel_bases


@dataclass
class NmRecoveryResult:
    """Result of NM tag recovery from converted-genome trials."""
    nm: int            # recovered NM value (-1 if unrecoverable)
    method: str        # how it was recovered
    trial_idx: int     # which trial matched (0-3), -1 if N/A
    all_trials: list   # all 4 trial distances
    confidence: str    # "high", "medium", "low"


def recover_bhmem_nm_from_trials(
    read,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> NmRecoveryResult:
    """Best-effort recovery of Bhmem ``NM:i`` from converted-genome CIGAR walks.

    No ``MD:Z`` used — purely CIGAR + bisulfite-converted FASTAs (CT/GA).

    **Strategy cascade** (each step covers reads the previous step missed):

    1. **Unique trial match:** if exactly one of the 4 trial distances equals ``NM:i``,
       that trial produced the alignment. High confidence (~92% R2, ~99% R1).
    2. **Multi trial match + Bhmem fold:** if multiple trials match ``NM:i``, pick among
       them using Bhmem single-record ordering (MAPQ/AS/NM/CIGAR-M — but on a fixed BAM
       these are tied, so any matching trial gives the correct distance). Medium confidence.
    3. **NM:i not in trials, NM=0:** heavily soft-clipped reads where BWA scored only the
       short aligned portion perfectly. Bhmem's ``modifySeqByCigar`` + PBAT strand-flip can
       corrupt stored SEQ for chimeric reads (multiple BWA SAM lines share a mutable
       ``originalSeq``), making trial recomputation wrong. We trust the tag. Medium confidence.
    4. **NM:i not in trials, NM>0:** same SEQ corruption, but NM>0. Rare (~0.3% of R2).
       Trust the tag. Low confidence.

    When ``NM:i`` is absent, falls back to ``min(trials)`` (the Bhmem fold result on a
    fixed BAM), tagged as ``"min_no_tag"`` / low confidence.

    Returns :class:`NmRecoveryResult`.
    """
    if read.is_unmapped or read.query_sequence is None:
        return NmRecoveryResult(0, "unmapped", -1, [], "high")

    trials = pbat_converted_genome_trial_distances(read, fasta_ct, fasta_ga)
    has_nm = read.has_tag("NM")
    nm_tag = int(read.get_tag("NM")) if has_nm else None

    if not trials:
        if has_nm:
            return NmRecoveryResult(nm_tag, "no_trials_trust_tag", -1, [], "low")
        return NmRecoveryResult(-1, "no_trials_no_tag", -1, [], "low")

    # No NM tag — fall back to min(trials)
    if not has_nm:
        mn = min(trials)
        idx = trials.index(mn)
        return NmRecoveryResult(mn, "min_no_tag", idx, trials, "low")

    # Find which trials match the NM tag
    matching = [i for i, d in enumerate(trials) if d == nm_tag]

    if len(matching) == 1:
        return NmRecoveryResult(nm_tag, "unique_trial", matching[0], trials, "high")

    if len(matching) > 1:
        return NmRecoveryResult(nm_tag, "multi_trial", matching[0], trials, "medium")

    # NM not in any trial — SEQ likely corrupted by modifySeqByCigar bug
    if nm_tag == 0:
        return NmRecoveryResult(0, "nm0_trust_tag", -1, trials, "medium")

    return NmRecoveryResult(nm_tag, "nm_nonzero_trust_tag", -1, trials, "low")


def recompute_bhmem_nm(
    read,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> int:
    """Convenience wrapper: return recovered NM value (int), or -1 if unrecoverable.

    Uses :func:`recover_bhmem_nm_from_trials`. No ``MD:Z``.
    """
    return recover_bhmem_nm_from_trials(read, fasta_ct, fasta_ga).nm


def main():
    ap = argparse.ArgumentParser(description="Bisulfite-corrected mismatch counts")
    ap.add_argument("reference_fasta", help="Genomic FASTA (e.g. mm10.fa, indexed .fai)")
    ap.add_argument("bam", help="BAM file")
    ap.add_argument("--max-reads", type=int, default=0, help="Limit reads (0 = all)")
    ap.add_argument(
        "--method",
        choices=("aligned_pairs", "cigar_mx"),
        default="aligned_pairs",
        help="aligned_pairs: get_aligned_pairs(matches_only=True). "
        "cigar_mx: only M/=/X columns (good for heavy soft-clip / indel BAMs like bhmem).",
    )
    args = ap.parse_args()

    counter = (
        count_mismatches_corrected_cigar_mx
        if args.method == "cigar_mx"
        else count_mismatches_corrected_aligned_pairs
    )

    fasta = pysam.FastaFile(args.reference_fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    n = 0
    tot_raw = tot_corr = tot_bs = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        raw, corr, bs = counter(read, fasta)
        tot_raw += raw
        tot_corr += corr
        tot_bs += bs
        n += 1
        if args.max_reads and n >= args.max_reads:
            break

    bam.close()
    fasta.close()

    print(f"method\t{args.method}")
    print(f"reads_used\t{n}")
    print(f"sum_raw_mismatch_positions\t{tot_raw}")
    print(f"sum_bisulfite_ignored_CT_or_GA\t{tot_bs}")
    print(f"sum_corrected_mismatch_positions\t{tot_corr}")
    print(
        "# Compare to NM:i tag: sum of per-read raw_mm should track NM for Bismark-like BAMs"
    )


if __name__ == "__main__":
    main()
