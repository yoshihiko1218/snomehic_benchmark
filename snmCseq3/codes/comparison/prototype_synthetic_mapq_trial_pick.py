#!/usr/bin/env python3
"""
Prototype: re-align reads with ``bwa mem -a`` against CT and GA indices to get
per-trial MAPQ at the BAM's fixed position.

For each read:
1. Convert query (R1: G->A, R2: C->T for -pbat)
2. Align converted query against full CT index and full GA index (``bwa mem -a``)
3. From BWA output, find all hits and their AS/MAPQ
4. For the BAM's fixed position, compute AS against each converted ref
5. Derive synthetic MAPQ: based on how the fixed-position AS compares to competing hits
6. Apply Bhmem cascade (MAPQ -> AS -> NM -> CIGAR-M) to pick trial
7. Compare picked trial's NM to the BAM's NM:i tag

Usage::

  python prototype_synthetic_mapq_trial_pick.py \\
    /path/to/bhmem.bam \\
    /path/to/Bisulfite_Genome \\
    --bwa /path/to/bwa \\
    --max-reads 1000
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile

import pysam
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import (
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
    _reverse_complement_dna,
)

DEFAULT_BISULFITE = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/Bisulfite_Genome"
DEFAULT_BWA = "/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/bwa"


def _convert_query_pbat(seq, is_read2):
    """PBAT conversion: R1 G->A, R2 C->T on the original FASTQ sequence."""
    if is_read2:
        return seq.upper().replace("C", "T")
    else:
        return seq.upper().replace("G", "A")


def _bwa_score_at_position(query_conv, ref_fasta, ref_contig, pos, cigar_tuples,
                           match=1, mismatch=4, gap_open=6, gap_extend=1):
    """Compute BWA-style alignment score at a fixed position with fixed CIGAR.

    Returns (AS, NM, n_matches, n_mismatches).
    """
    if ref_contig not in ref_fasta.references:
        return None

    q = 0
    r = pos
    n_match = 0
    n_mismatch = 0
    gap_penalty = 0

    for op, length in cigar_tuples:
        if op in (0, 7, 8):  # M, =, X
            ref_span = ref_fasta.fetch(ref_contig, r, r + length).upper()
            for i in range(length):
                if q < len(query_conv) and i < len(ref_span):
                    if query_conv[q] == ref_span[i]:
                        n_match += 1
                    else:
                        n_mismatch += 1
                q += 1
                r += 1
        elif op == 1:  # I
            gap_penalty += gap_open + gap_extend * length
            q += length
        elif op in (2, 3):  # D, N
            gap_penalty += gap_open + gap_extend * length
            r += length
        elif op == 4:  # S
            q += length
        elif op == 5:  # H
            pass

    AS = n_match * match - n_mismatch * mismatch - gap_penalty
    NM = n_mismatch + sum(l for o, l in cigar_tuples if o == 1) + sum(l for o, l in cigar_tuples if o in (2, 3))
    return AS, NM, n_match, n_mismatch


def _run_bwa_mem_all_hits(bwa_path, index_fa, query_seq, query_name="read"):
    """Run bwa mem -a on a single read, return list of (chrom, pos, mapq, AS, cigar, flag)."""
    hits = []
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fq", delete=False) as fq:
        fq.write(f"@{query_name}\n{query_seq}\n+\n{'I' * len(query_seq)}\n")
        fq_path = fq.name

    try:
        result = subprocess.run(
            [bwa_path, "mem", "-a", "-t", "1", index_fa, fq_path],
            capture_output=True, text=True, timeout=30,
        )
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4:  # unmapped
                continue
            chrom = fields[2]
            # Strip conversion suffix from contig name
            chrom_clean = chrom.replace("_CT_converted", "").replace("_GA_converted", "")
            pos = int(fields[3]) - 1  # SAM is 1-based, convert to 0-based
            mapq = int(fields[4])
            cigar = fields[5]

            # Parse tags for AS
            as_val = None
            for f in fields[11:]:
                if f.startswith("AS:i:"):
                    as_val = int(f.split(":")[2])
                    break

            hits.append({
                "chrom": chrom_clean,
                "chrom_raw": chrom,
                "pos": pos,
                "mapq": mapq,
                "AS": as_val,
                "cigar": cigar,
                "flag": flag,
            })
    finally:
        os.unlink(fq_path)

    return hits


def _estimate_mapq_at_position(target_AS, all_hits_AS, n_hits):
    """Rough BWA-MEM-style MAPQ estimate.

    BWA-MEM MAPQ ≈ -10*log10(P(wrong)), where P(wrong) depends on:
    - gap between best AS and second-best AS
    - number of equally-good hits

    Simplified model:
    - If target_AS is strictly best and gap >= 20: MAPQ = 60
    - If target_AS == best and unique: scale by gap to second
    - If target_AS < best: MAPQ = 0
    - If multiple hits at same AS as target: MAPQ = 0-3
    """
    if not all_hits_AS or target_AS is None:
        return 0

    best_AS = max(all_hits_AS)

    if target_AS < best_AS:
        # Fixed position is not the best hit
        return 0

    # target_AS == best_AS
    n_at_best = sum(1 for a in all_hits_AS if a == target_AS)
    if n_at_best > 1:
        # Multiple hits at same score
        return min(3, max(0, int(-10 * np.log10(1 - 1.0 / n_at_best + 1e-10))))

    # Unique best hit — MAPQ depends on gap to second best
    sorted_as = sorted(all_hits_AS, reverse=True)
    if len(sorted_as) < 2:
        return 60  # no second hit

    gap = sorted_as[0] - sorted_as[1]
    if gap >= 30:
        return 60
    if gap >= 20:
        return 50
    if gap >= 10:
        return 40
    if gap >= 5:
        return 25
    if gap >= 1:
        return 10
    return 3


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", help="Bhmem BAM file")
    ap.add_argument("bisulfite_genome", nargs="?", default=DEFAULT_BISULFITE,
                    help=f"Bisulfite_Genome dir (default: {DEFAULT_BISULFITE})")
    ap.add_argument("--bwa", default=DEFAULT_BWA, help="Path to bwa")
    ap.add_argument("--max-reads", type=int, default=1000)
    ap.add_argument("-o", "--output-prefix", default="synthetic_mapq_test")
    args = ap.parse_args()

    ct_fa_path = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_fa_path = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")

    # Check BWA indices exist
    for fa in (ct_fa_path, ga_fa_path):
        if not os.path.isfile(fa + ".bwt"):
            print(f"ERROR: BWA index missing for {fa} — run 'bwa index {fa}' first", file=sys.stderr)
            sys.exit(1)

    ct_fa = pysam.FastaFile(ct_fa_path)
    ga_fa = pysam.FastaFile(ga_fa_path)
    bam = pysam.AlignmentFile(args.bam, "rb")

    n = 0
    results = []

    # Stats
    correct_synth = 0
    correct_min = 0
    total_with_nm = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None or not read.has_tag("NM"):
            continue

        n += 1
        if n > args.max_reads:
            break
        if n % 100 == 0:
            print(f"  {n} reads...", flush=True)

        nm_tag = int(read.get_tag("NM"))
        is_r2 = bool(read.is_paired and read.is_read2)
        bam_chrom = read.reference_name
        bam_pos = read.reference_start
        bam_cigar = read.cigartuples
        bam_mapq = read.mapping_quality

        # Convert query as Bhmem would
        seq = read.query_sequence.upper()
        query_conv = _convert_query_pbat(seq, is_r2)

        # Also try RC conversion (for reverse-strand stored reads)
        query_conv_v2 = _reverse_complement_dna(_convert_query_pbat(_reverse_complement_dna(seq), is_r2))

        # Run bwa mem -a against both indices
        ct_hits = _run_bwa_mem_all_hits(args.bwa, ct_fa_path, query_conv, read.query_name)
        ga_hits = _run_bwa_mem_all_hits(args.bwa, ga_fa_path, query_conv, read.query_name)

        # Also try v2 conversion
        ct_hits_v2 = _run_bwa_mem_all_hits(args.bwa, ct_fa_path, query_conv_v2, read.query_name)
        ga_hits_v2 = _run_bwa_mem_all_hits(args.bwa, ga_fa_path, query_conv_v2, read.query_name)

        # For each trial, compute AS at the BAM's fixed position and synthetic MAPQ
        trials = []
        for conv_label, conv_q, fa, fa_path, hits_list in [
            ("CT_v1", query_conv, ct_fa, ct_fa_path,  ct_hits),
            ("GA_v1", query_conv, ga_fa, ga_fa_path,  ga_hits),
            ("CT_v2", query_conv_v2, ct_fa, ct_fa_path, ct_hits_v2),
            ("GA_v2", query_conv_v2, ga_fa, ga_fa_path, ga_hits_v2),
        ]:
            genome = conv_label[:2]
            c_name = bisulfite_converted_contig_name(fa, bam_chrom, genome)
            if c_name is None:
                continue

            score_result = _bwa_score_at_position(
                conv_q, fa, c_name, bam_pos, bam_cigar
            )
            if score_result is None:
                continue

            as_fixed, nm_fixed, n_match, n_mismatch = score_result

            # Collect all AS from bwa hits for MAPQ estimation
            all_as = [h["AS"] for h in hits_list if h["AS"] is not None]

            synth_mapq = _estimate_mapq_at_position(as_fixed, all_as, len(hits_list))

            trials.append({
                "label": conv_label,
                "AS": as_fixed,
                "NM": nm_fixed,
                "synth_mapq": synth_mapq,
                "n_bwa_hits": len(hits_list),
            })

        if not trials:
            continue

        total_with_nm += 1

        # Pick by Bhmem cascade: MAPQ -> AS -> NM -> (skip CIGAR-M, same for all)
        def trial_key(t):
            return (t["synth_mapq"], t["AS"], -t["NM"])

        best_synth = max(trials, key=trial_key)
        best_min = min(trials, key=lambda t: t["NM"])

        if best_synth["NM"] == nm_tag:
            correct_synth += 1
        if best_min["NM"] == nm_tag:
            correct_min += 1

        results.append({
            "name": read.query_name,
            "is_r2": is_r2,
            "nm_tag": nm_tag,
            "bam_mapq": bam_mapq,
            "synth_pick_nm": best_synth["NM"],
            "synth_pick_mapq": best_synth["synth_mapq"],
            "synth_pick_label": best_synth["label"],
            "min_pick_nm": best_min["NM"],
            "all_trials": [(t["label"], t["synth_mapq"], t["AS"], t["NM"]) for t in trials],
        })

    bam.close()
    ct_fa.close()
    ga_fa.close()

    # Report
    print(f"\n{'='*60}")
    print(f"Synthetic MAPQ trial pick prototype ({total_with_nm} reads with NM tag)")
    print(f"{'='*60}\n")

    if total_with_nm == 0:
        print("No reads processed.")
        return

    print(f"Synthetic MAPQ pick matches NM tag: {correct_synth}/{total_with_nm} ({100*correct_synth/total_with_nm:.1f}%)")
    print(f"min(trials) matches NM tag:         {correct_min}/{total_with_nm} ({100*correct_min/total_with_nm:.1f}%)")

    # Break down by R1/R2
    for mate_label, mate_flag in [("R1", False), ("R2", True)]:
        subset = [r for r in results if r["is_r2"] == mate_flag]
        if not subset:
            continue
        ns = len(subset)
        cs = sum(1 for r in subset if r["synth_pick_nm"] == r["nm_tag"])
        cm = sum(1 for r in subset if r["min_pick_nm"] == r["nm_tag"])
        print(f"\n  {mate_label} ({ns} reads):")
        print(f"    Synthetic MAPQ: {cs}/{ns} ({100*cs/ns:.1f}%)")
        print(f"    min(trials):    {cm}/{ns} ({100*cm/ns:.1f}%)")

    # Write TSV
    tsv_path = f"{args.output_prefix}.per_read.tsv"
    with open(tsv_path, "w") as f:
        f.write("read_name\tis_r2\tnm_tag\tbam_mapq\tsynth_pick_nm\tsynth_pick_mapq\t"
                "synth_pick_label\tmin_pick_nm\tsynth_correct\tmin_correct\tall_trials\n")
        for r in results:
            sc = 1 if r["synth_pick_nm"] == r["nm_tag"] else 0
            mc = 1 if r["min_pick_nm"] == r["nm_tag"] else 0
            trials_str = ";".join(f"{t[0]}:MQ={t[1]},AS={t[2]},NM={t[3]}" for t in r["all_trials"])
            f.write(f"{r['name']}\t{int(r['is_r2'])}\t{r['nm_tag']}\t{r['bam_mapq']}\t"
                    f"{r['synth_pick_nm']}\t{r['synth_pick_mapq']}\t{r['synth_pick_label']}\t"
                    f"{r['min_pick_nm']}\t{sc}\t{mc}\t{trials_str}\n")
    print(f"\nWrote {tsv_path}")


if __name__ == "__main__":
    main()
