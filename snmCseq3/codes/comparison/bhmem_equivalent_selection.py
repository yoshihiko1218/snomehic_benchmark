#!/usr/bin/env python3
"""
Bhmem-equivalent **alignment record selection** and **NM** from CIGAR + bisulfite **converted** FASTA.

Mirrors ``edu.mit.compbio.bisulfitehic.mapping.Bhmem`` (bisulfitehic jar) for:

- ``comparingSamRecord(SAMRecord, SAMRecord)`` — single-end merge within one index stream (MAPQ,
  ``AS``, ``NM``, then CIGAR ``M`` span).
- ``comparingSamRecordPbat(Pair, Pair, regionsEnzyme)`` — paired-end choice for ``-pbat`` +
  ``-nonDirectional`` including **restriction enzyme** tie-break when ``-enzymeList`` is loaded
  (:class:`EnzymeRegionIndex`).

**CIGAR tie-break:** matches htsjdk ``CigarUtil.cigarArrayFromString``: count **per-base** ``M``
operators only (BAM op ``M`` / ``0``); ``=`` and ``X`` are **not** counted, same as Java
``if (cigar == 'M')``.

**NM from converted reference:** uses :func:`bisulfite_corrected_mismatch.count_nm_style_edit_distance_converted_explicit`
with PBAT converted-query variants (R1 G→A, R2 C→T and the two orientations from
``_pbat_converted_query_variants``). No ``MD:Z``.

**When ``AS`` is unavailable** (e.g. some Yap/Bismark records): each trial uses
``AS_proxy = -NM_trial`` so that higher proxy AS favors lower recomputed edit distance, consistent
with the ordering after MAPQ when true BWA ``AS`` would correlate with alignment quality. For an
exact replay of Bhmem, supply real per-trial ``AS`` from jbwa (same scoring as BWA).

**Trial disambiguation** when several conversions tie on a fixed CIGAR: see
:func:`pick_pbat_single_trial` (``margin_at_min_dist``, NM-guided modes) and
:func:`pick_pbat_single_trial_margin_at_min_dist`.

Reference: ``Bhmem.java`` lines 717--736 (single), 817--864 (paired PBAT with ``-enzymeList``).
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Callable, Iterable, TypeVar

import pysam

from bisulfite_corrected_mismatch import (
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
    _pbat_converted_query_variants,
)


# Bhmem expands enzyme overlap query by ±50 on **1-based inclusive** alignment ends (SAM / htsjdk).
_ENZYME_PAD = 50


@dataclass
class EnzymeRegionIndex:
    """Restriction intervals per reference name (0-based half-open), for Bhmem-style overlap tests."""

    _starts: dict[str, list[int]] = field(default_factory=dict)
    _ends: dict[str, list[int]] = field(default_factory=dict)

    @classmethod
    def from_bed_file(cls, path: str) -> EnzymeRegionIndex:
        """Load BED / bedgraph (chr, start, end tab-separated). Coordinates 0-based half-open like BED."""
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 3:
                    continue
                chrom, a, b = p[0], p[1], p[2]
                try:
                    s, e = int(a), int(b)
                except ValueError:
                    continue
                if e <= s:
                    continue
                by_chrom.setdefault(chrom, []).append((s, e))
        idx = cls()
        for chrom, ivs in by_chrom.items():
            ivs.sort()
            merged: list[tuple[int, int]] = []
            for s, e in ivs:
                if not merged or s > merged[-1][1]:
                    merged.append((s, e))
                else:
                    pl, pr = merged[-1]
                    merged[-1] = (pl, max(pr, e))
            idx._starts[chrom] = [x[0] for x in merged]
            idx._ends[chrom] = [x[1] for x in merged]
        return idx

    def has_chrom(self, chrom: str) -> bool:
        return chrom in self._starts

    def is_empty(self) -> bool:
        return not self._starts

    @property
    def n_chromosomes(self) -> int:
        return len(self._starts)

    def overlaps_padded_alignment(self, read: pysam.AlignedSegment) -> bool:
        """Match Bhmem: ``overlappers(alignmentStart-50, alignmentEnd+50)`` (htsjdk 1-based coords)."""
        if read.is_unmapped:
            return False
        chrom = read.reference_name
        if chrom not in self._starts:
            return False
        s0 = int(read.reference_start)
        e0_excl = int(read.reference_end)
        # 1-based inclusive SAM coords: start = s0+1, end = e0_excl (pysam exclusive end == 1-based last + 1)
        j_lo_1 = s0 + 1 - _ENZYME_PAD
        j_hi_1 = e0_excl + _ENZYME_PAD
        q0 = max(0, j_lo_1 - 1)
        q1 = max(q0, j_hi_1)
        return self._overlaps_half_open(chrom, q0, q1)

    def _overlaps_half_open(self, chrom: str, q0: int, q1: int) -> bool:
        ss = self._starts.get(chrom)
        ee = self._ends.get(chrom)
        if not ss:
            return False
        i = bisect.bisect_right(ss, q1) - 1
        while i >= 0:
            if ee[i] > q0:
                return True
            if ss[i] < q0:
                break
            i -= 1
        return False


def _safe_int_tag(read: pysam.AlignedSegment, tag: str) -> int | None:
    if not read.has_tag(tag):
        return None
    v = read.get_tag(tag)
    if v is None:
        return None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def cigar_m_len_bhmem(read: pysam.AlignedSegment) -> int:
    """Per-base count of ``M`` in htsjdk-expanded CIGAR sense (not ``=`` / ``X``)."""
    ct = read.cigartuples
    if not ct:
        return 0
    return sum(length for op, length in ct if op == 0)


def bhmem_prefer_second_single(
    mq1: int,
    as1: int,
    nm1: int,
    cm1: int,
    mq2: int,
    as2: int,
    nm2: int,
    cm2: int,
) -> bool:
    """Return True if the second record should win (Bhmem ``comparingSamRecord``)."""
    if mq2 > mq1:
        return True
    if mq2 < mq1:
        return False
    if as2 > as1:
        return True
    if as2 < as1:
        return False
    if nm2 < nm1:
        return True
    if nm2 > nm1:
        return False
    return cm2 > cm1


_T = TypeVar("_T")


def fold_best_single(
    items: Iterable[_T],
    key: Callable[[_T], tuple[int, int, int, int]],
) -> _T | None:
    """Pick one item by pairwise Bhmem single-end ordering (left = current best, right = challenger)."""
    it = iter(items)
    try:
        best = next(it)
    except StopIteration:
        return None
    bmq, bas, bnm, bcm = key(best)
    for x in it:
        mq, as_, nm, cm = key(x)
        if bhmem_prefer_second_single(bmq, bas, bnm, bcm, mq, as_, nm, cm):
            best = x
            bmq, bas, bnm, bcm = mq, as_, nm, cm
    return best


@dataclass(frozen=True)
class PbatTrialLabel:
    genome: str  # "CT" or "GA"
    qtag: str  # "seq_conv" or "rc_conv_rc"


@dataclass(frozen=True)
class PbatSingleTrial:
    label: PbatTrialLabel
    dist: int
    fasta: pysam.FastaFile
    ref_contig: str
    query_converted: str


def enumerate_pbat_single_trials(
    read: pysam.AlignedSegment,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> list[PbatSingleTrial]:
    """All valid CT/GA × PBAT query-orientation trials with recomputed NM-style distance."""
    if read.is_unmapped or read.query_sequence is None:
        return []
    c_ct = bisulfite_converted_contig_name(fasta_ct, read.reference_name, "CT")
    c_ga = bisulfite_converted_contig_name(fasta_ga, read.reference_name, "GA")
    if c_ct is None or c_ga is None:
        return []

    out: list[PbatSingleTrial] = []
    for conv, fa, cname in (
        ("CT", fasta_ct, c_ct),
        ("GA", fasta_ga, c_ga),
    ):
        for qtag, qconv in _pbat_converted_query_variants(read):
            d = count_nm_style_edit_distance_converted_explicit(
                read, fa, ref_contig=cname, query_converted=qconv
            )
            if d < 0:
                continue
            out.append(
                PbatSingleTrial(
                    label=PbatTrialLabel(conv, qtag),
                    dist=d,
                    fasta=fa,
                    ref_contig=cname,
                    query_converted=qconv,
                )
            )
    return out


def recompute_nm_bhmem_style_single_pbat(
    read: pysam.AlignedSegment,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
    *,
    as_per_trial: dict[tuple[str, str], int] | None = None,
) -> tuple[int, PbatTrialLabel | None]:
    """Pick a trial using the **same ordering as** Bhmem ``comparingSamRecord`` (MAPQ→AS→NM→CIGAR-M).

    For each trial we build the tuple ``(MAPQ, AS, NM, cigar_M_len)`` and fold with
    :func:`bhmem_prefer_second_single` (higher MQ/AS, lower NM, more ``M`` bases win).

    **Fields on a fixed BAM record**

    - **MAPQ** — from the read; **identical** for all four trials (one alignment).
    - **AS** — if ``as_per_trial`` is set, **per-trial** BWA scores (jbwa). Else ``AS:i`` from the
      BAM if present: **identical** across trials. Else **``-dist``** per trial so “higher AS”
      matches lower recomputed edit distance (Bhmem would have a different ``AS`` per BWA line).
    - **NM** — we use **recomputed** ``dist`` for that trial **not** ``NM:i`` on the read for every
      trial: the tag is the **winner’s** single value; repeating it would make NM tie on all
      trials and **remove** the discriminator Bhmem has when comparing **different** SAM lines.
    - **CIGAR M span** — :func:`cigar_m_len_bhmem`; **identical** across trials (same CIGAR).

    So nothing is “skipped”: MQ/AS/M are applied, but only **NM (as recomputed per trial)** and
    optionally **per-trial AS** actually separate trials without jbwa.

    ``as_per_trial`` maps ``(genome, qtag)`` e.g. ``("CT", "seq_conv")`` to BWA ``AS`` for that
    realignment.
    """
    trials = enumerate_pbat_single_trials(read, fasta_ct, fasta_ga)
    if not trials:
        return -1, None

    mq = int(read.mapping_quality)
    cm = cigar_m_len_bhmem(read)
    as_common = _safe_int_tag(read, "AS")

    def metrics(tr: PbatSingleTrial) -> tuple[int, int, int, int]:
        # Order matches Bhmem comparingSamRecord: (MQ, AS, NM, cigarM); fold uses bhmem_prefer_second_single.
        key = (tr.label.genome, tr.label.qtag)
        if as_per_trial is not None and key in as_per_trial:
            as_ = int(as_per_trial[key])
        elif as_common is not None:
            as_ = as_common
        else:
            as_ = -tr.dist
        return mq, as_, tr.dist, cm

    best = fold_best_single(trials, metrics)
    assert best is not None
    return best.dist, best.label


def _tiebreak_pbat_trials_bhmem(
    trials: list[PbatSingleTrial],
    read: pysam.AlignedSegment,
    *,
    as_per_trial: dict[tuple[str, str], int] | None = None,
) -> PbatSingleTrial:
    """Bhmem single-record ordering among a **subset** of trials (MAPQ, AS, NM, CIGAR-M)."""
    mq = int(read.mapping_quality)
    cm = cigar_m_len_bhmem(read)
    as_common = _safe_int_tag(read, "AS")

    def metrics(tr: PbatSingleTrial) -> tuple[int, int, int, int]:
        key = (tr.label.genome, tr.label.qtag)
        if as_per_trial is not None and key in as_per_trial:
            as_ = int(as_per_trial[key])
        elif as_common is not None:
            as_ = as_common
        else:
            as_ = -tr.dist
        return mq, as_, tr.dist, cm

    best = fold_best_single(trials, metrics)
    assert best is not None
    return best


def pick_pbat_single_trial_margin_at_min_dist(
    read: pysam.AlignedSegment,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
    *,
    as_per_trial: dict[tuple[str, str], int] | None = None,
) -> tuple[int, PbatTrialLabel | None]:
    """Among trials with **global** minimum recomputed distance, pick the one with largest **margin**.

    Margin for trial ``t`` is ``min(dist(u) for u != t) - dist(t)``. Intuition: the true conversion
    often separates clearly from the other three; when two trials tie at the minimum with a tight
    runner-up, Bhmem-style tie-break applies.

    This uses **no** ``NM:i`` — suitable when the tag is missing or not comparable (e.g. some Yap).
    """
    trials = enumerate_pbat_single_trials(read, fasta_ct, fasta_ga)
    if not trials:
        return -1, None
    m = min(t.dist for t in trials)
    pool = [t for t in trials if t.dist == m]
    if len(pool) == 1:
        return pool[0].dist, pool[0].label

    def margin(tr: PbatSingleTrial) -> int:
        ods = [u.dist for u in trials if u is not tr]
        return (min(ods) - tr.dist) if ods else 10**9

    max_mar = max(margin(t) for t in pool)
    pool2 = [t for t in pool if margin(t) == max_mar]
    best = _tiebreak_pbat_trials_bhmem(pool2, read, as_per_trial=as_per_trial)
    return best.dist, best.label


def pick_pbat_single_trial(
    read: pysam.AlignedSegment,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
    *,
    strategy: str = "bhmem_fold",
    as_per_trial: dict[tuple[str, str], int] | None = None,
    nm_hint: int | None = None,
) -> tuple[int, PbatTrialLabel | None, str]:
    """Choose one CT/GA × PBAT-orientation trial by **strategy** (returns distance, label, reason).

    Strategies:

    - ``bhmem_fold`` — :func:`recompute_nm_bhmem_style_single_pbat` (MAPQ→AS→NM→M; on a fixed BAM
      usually **min distance** plus CIGAR-M tie-break).
    - ``margin_at_min_dist`` — :func:`pick_pbat_single_trial_margin_at_min_dist`.
    - ``unique_nm_match_else_bhmem`` — if ``NM`` (``nm_hint`` or ``NM:i``) matches **exactly one**
      trial’s distance, use that trial; else ``bhmem_fold``. **Bhmem / BWA calibration only** —
      circular if you use ``NM:i`` to define “truth” and the same tag to pick.
    - ``unique_nm_match_else_nearest`` — unique ``NM`` hit else minimize ``(abs(d-NM), d)`` among
      trials, Bhmem tie-break on ties.
    - ``nearest_nm`` — minimize ``(abs(d-NM), d)``; requires ``nm_hint`` or ``NM:i``; else falls
      back to ``bhmem_fold`` with reason ``no_nm_tag``.
    """
    nm = nm_hint if nm_hint is not None else _safe_int_tag(read, "NM")

    if strategy == "bhmem_fold":
        d, lab = recompute_nm_bhmem_style_single_pbat(
            read, fasta_ct, fasta_ga, as_per_trial=as_per_trial
        )
        return d, lab, "bhmem_fold"

    if strategy == "margin_at_min_dist":
        d, lab = pick_pbat_single_trial_margin_at_min_dist(
            read, fasta_ct, fasta_ga, as_per_trial=as_per_trial
        )
        return d, lab, "margin_at_min_dist"

    trials = enumerate_pbat_single_trials(read, fasta_ct, fasta_ga)
    if not trials:
        return -1, None, "no_trials"

    if strategy == "unique_nm_match_else_bhmem":
        if nm is None:
            d, lab = recompute_nm_bhmem_style_single_pbat(
                read, fasta_ct, fasta_ga, as_per_trial=as_per_trial
            )
            return d, lab, "fallback_bhmem_no_nm_tag"
        exact = [t for t in trials if t.dist == nm]
        if len(exact) == 1:
            return exact[0].dist, exact[0].label, "unique_nm_match"
        d, lab = recompute_nm_bhmem_style_single_pbat(
            read, fasta_ct, fasta_ga, as_per_trial=as_per_trial
        )
        return d, lab, "fallback_bhmem"

    if strategy == "nearest_nm":
        if nm is None:
            d, lab = recompute_nm_bhmem_style_single_pbat(
                read, fasta_ct, fasta_ga, as_per_trial=as_per_trial
            )
            return d, lab, "no_nm_tag"
        best_key = min((abs(t.dist - nm), t.dist) for t in trials)
        pool = [t for t in trials if (abs(t.dist - nm), t.dist) == best_key]
        best = _tiebreak_pbat_trials_bhmem(pool, read, as_per_trial=as_per_trial)
        return best.dist, best.label, "nearest_nm"

    if strategy == "unique_nm_match_else_nearest":
        if nm is None:
            d, lab = recompute_nm_bhmem_style_single_pbat(
                read, fasta_ct, fasta_ga, as_per_trial=as_per_trial
            )
            return d, lab, "fallback_bhmem_no_nm_tag"
        exact = [t for t in trials if t.dist == nm]
        if len(exact) == 1:
            return exact[0].dist, exact[0].label, "unique_nm_match"
        best_key = min((abs(t.dist - nm), t.dist) for t in trials)
        pool = [t for t in trials if (abs(t.dist - nm), t.dist) == best_key]
        best = _tiebreak_pbat_trials_bhmem(pool, read, as_per_trial=as_per_trial)
        return best.dist, best.label, "nearest_nm"

    raise ValueError(f"unknown strategy: {strategy!r}")


@dataclass(frozen=True)
class PbatPairAssignment:
    conv1: str
    conv2: str
    qtag1: str
    qtag2: str
    d1: int
    d2: int


def _same_contig(a: pysam.AlignedSegment, b: pysam.AlignedSegment) -> bool:
    return str(a.reference_name).lower() == str(b.reference_name).lower()


def _enzyme_branch_prefers_candidate_b(
    enzyme: EnzymeRegionIndex | None,
    ra1: pysam.AlignedSegment,
    ra2: pysam.AlignedSegment,
    rb1: pysam.AlignedSegment,
    rb2: pysam.AlignedSegment,
) -> bool:
    """Bhmem ``returnR2`` from restriction-enzyme logic (lines 749--763)."""
    if enzyme is None or enzyme.is_empty():
        return False
    if not (
        enzyme.has_chrom(rb1.reference_name) and enzyme.has_chrom(rb2.reference_name)
    ):
        return False
    if not enzyme.has_chrom(ra1.reference_name) or not enzyme.has_chrom(
        ra2.reference_name
    ):
        return True
    b_hit = enzyme.overlaps_padded_alignment(rb1) and enzyme.overlaps_padded_alignment(
        rb2
    )
    if not b_hit:
        return False
    a_hit1 = enzyme.overlaps_padded_alignment(ra1)
    a_hit2 = enzyme.overlaps_padded_alignment(ra2)
    return (not a_hit1) or (not a_hit2)


def bhmem_prefer_second_pair_pbat(
    mq_a1: int,
    mq_a2: int,
    as_a1: int,
    as_a2: int,
    nm_a1: int,
    nm_a2: int,
    cm_a1: int,
    cm_a2: int,
    chr_same_a: bool,
    mq_b1: int,
    mq_b2: int,
    as_b1: int,
    as_b2: int,
    nm_b1: int,
    nm_b2: int,
    cm_b1: int,
    cm_b2: int,
    chr_same_b: bool,
    enzyme: EnzymeRegionIndex | None,
    ra1: pysam.AlignedSegment,
    ra2: pysam.AlignedSegment,
    rb1: pysam.AlignedSegment,
    rb2: pysam.AlignedSegment,
) -> bool:
    """Paired-end Bhmem ``comparingSamRecordPbat(Pair, Pair, regionsEnzyme)`` (lines 817--864)."""
    if (mq_b1 > 0 and mq_b2 > 0) and (mq_a1 == 0 or mq_a2 == 0):
        return True

    sa = mq_a1 + mq_a2
    sb = mq_b1 + mq_b2
    if sb > sa:
        return True
    if sb < sa:
        return False
    if (not chr_same_a) and chr_same_b:
        return True
    if chr_same_a and (not chr_same_b):
        return False

    if _enzyme_branch_prefers_candidate_b(enzyme, ra1, ra2, rb1, rb2):
        return True

    sas = as_a1 + as_a2
    sbs = as_b1 + as_b2
    if sbs > sas:
        return True
    if sbs < sas:
        return False
    nma = nm_a1 + nm_a2
    nmb = nm_b1 + nm_b2
    if nmb < nma:
        return True
    if nmb > nma:
        return False
    cma = cm_a1 + cm_a2
    cmb = cm_b1 + cm_b2
    return cmb > cma


def fold_best_pair(
    assignments: list[PbatPairAssignment],
    read1: pysam.AlignedSegment,
    read2: pysam.AlignedSegment,
    *,
    enzyme: EnzymeRegionIndex | None = None,
) -> PbatPairAssignment | None:
    """Fold assignments with Bhmem ``comparingSamRecordPbat`` (optional ``-enzymeList`` index).

    Uses per-assignment ``NM`` = recomputed converted-genome distance and ``AS`` = tag when present
    else ``-d`` per mate. Mate coordinates are identical across assignments, so the enzyme block
    rarely changes ordering when recomputing from a fixed BAM; it matters when candidate pairs
    differ in locus (e.g. full jbwa replay).
    """
    if not assignments:
        return None
    mq1 = int(read1.mapping_quality)
    mq2 = int(read2.mapping_quality)
    cm1 = cigar_m_len_bhmem(read1)
    cm2 = cigar_m_len_bhmem(read2)
    as1_tag = _safe_int_tag(read1, "AS")
    as2_tag = _safe_int_tag(read2, "AS")

    def vec(ass: PbatPairAssignment):
        as1 = as1_tag if as1_tag is not None else -ass.d1
        as2 = as2_tag if as2_tag is not None else -ass.d2
        return (
            mq1,
            mq2,
            as1,
            as2,
            ass.d1,
            ass.d2,
            cm1,
            cm2,
            _same_contig(read1, read2),
        )

    best = assignments[0]
    va = vec(best)
    for cand in assignments[1:]:
        vb = vec(cand)
        if bhmem_prefer_second_pair_pbat(
            va[0],
            va[1],
            va[2],
            va[3],
            va[4],
            va[5],
            va[6],
            va[7],
            va[8],
            vb[0],
            vb[1],
            vb[2],
            vb[3],
            vb[4],
            vb[5],
            vb[6],
            vb[7],
            vb[8],
            enzyme,
            read1,
            read2,
            read1,
            read2,
        ):
            best = cand
            va = vb
    return best


def enumerate_pbat_pair_assignments(
    read1: pysam.AlignedSegment,
    read2: pysam.AlignedSegment,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
) -> list[PbatPairAssignment]:
    pair_index = (("CT", "CT"), ("GA", "GA"), ("CT", "GA"), ("GA", "CT"))
    out: list[PbatPairAssignment] = []
    for idx1, idx2 in pair_index:
        fa1 = fasta_ct if idx1 == "CT" else fasta_ga
        fa2 = fasta_ct if idx2 == "CT" else fasta_ga
        c1 = bisulfite_converted_contig_name(fa1, read1.reference_name, idx1)
        c2 = bisulfite_converted_contig_name(fa2, read2.reference_name, idx2)
        if c1 is None or c2 is None:
            continue
        for t1, q1 in _pbat_converted_query_variants(read1):
            d1 = count_nm_style_edit_distance_converted_explicit(
                read1, fa1, ref_contig=c1, query_converted=q1
            )
            if d1 < 0:
                continue
            for t2, q2 in _pbat_converted_query_variants(read2):
                d2 = count_nm_style_edit_distance_converted_explicit(
                    read2, fa2, ref_contig=c2, query_converted=q2
                )
                if d2 < 0:
                    continue
                out.append(
                    PbatPairAssignment(idx1, idx2, t1, t2, d1, d2)
                )
    return out


def recompute_nm_bhmem_style_pair_pbat_nd(
    read1: pysam.AlignedSegment,
    read2: pysam.AlignedSegment,
    fasta_ct: pysam.FastaFile,
    fasta_ga: pysam.FastaFile,
    *,
    enzyme: EnzymeRegionIndex | None = None,
) -> tuple[int, int, PbatPairAssignment | None]:
    """Non-directional PBAT pair: Bhmem-style fold (optional enzyme index); return (d1, d2, pick)."""
    if (
        read1.is_unmapped
        or read2.is_unmapped
        or read1.query_sequence is None
        or read2.query_sequence is None
    ):
        return -1, -1, None
    ass = enumerate_pbat_pair_assignments(read1, read2, fasta_ct, fasta_ga)
    if not ass:
        return -1, -1, None
    best = fold_best_pair(ass, read1, read2, enzyme=enzyme)
    assert best is not None
    return best.d1, best.d2, best
