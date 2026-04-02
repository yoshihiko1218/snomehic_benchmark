#!/usr/bin/env python3
"""
Batched synthetic MAPQ prototype: re-align reads with ``bwa mem -a`` against
CT and GA indices to get per-trial hit landscape, then pick trial using
Bhmem-style cascade with synthetic MAPQ.

Steps:
1. Read BAM, extract converted queries into 4 FASTQ files (2 conversions x 2 indices)
2. Run ``bwa mem -a`` once per FASTQ (4 total BWA runs, fast)
3. Parse all hits, group by read
4. For each read's BAM position, compute AS against each converted ref
5. Estimate MAPQ from hit landscape
6. Pick trial with Bhmem cascade, compare to NM:i tag
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bisulfite_corrected_mismatch import (
    bisulfite_converted_contig_name,
    _reverse_complement_dna,
)

DEFAULT_BISULFITE = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/Bisulfite_Genome"
DEFAULT_BWA = "/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/bwa"


def _convert_query_pbat(seq, is_read2):
    if is_read2:
        return seq.upper().replace("C", "T")
    else:
        return seq.upper().replace("G", "A")


def _bwa_score_at_pos(query_conv, ref_fasta, ref_contig, pos, cigar_tuples,
                      match=1, mismatch=4, gap_open=6, gap_extend=1):
    """Compute BWA AS and NM at a fixed position."""
    if ref_contig not in ref_fasta.references:
        return None
    q = 0
    r = pos
    n_match = 0
    n_mismatch = 0
    gap_pen = 0
    for op, length in cigar_tuples:
        if op in (0, 7, 8):
            ref_span = ref_fasta.fetch(ref_contig, r, r + length).upper()
            for i in range(length):
                if q < len(query_conv) and i < len(ref_span):
                    if query_conv[q] == ref_span[i]:
                        n_match += 1
                    else:
                        n_mismatch += 1
                q += 1
                r += 1
        elif op == 1:
            gap_pen += gap_open + gap_extend * length
            q += length
        elif op in (2, 3):
            gap_pen += gap_open + gap_extend * length
            r += length
        elif op == 4:
            q += length
        elif op == 5:
            pass
    AS = n_match * match - n_mismatch * mismatch - gap_pen
    NM = n_mismatch + sum(l for o, l in cigar_tuples if o == 1) + sum(l for o, l in cigar_tuples if o in (2, 3))
    return AS, NM


def _estimate_mapq(target_AS, all_AS_values):
    """Simplified BWA-MEM MAPQ from AS landscape."""
    if not all_AS_values or target_AS is None:
        return 0
    best = max(all_AS_values)
    if target_AS < best:
        return 0
    n_at_best = sum(1 for a in all_AS_values if a == target_AS)
    if n_at_best > 1:
        return min(3, max(0, int(-10 * np.log10(1.0 - 1.0 / n_at_best + 1e-10))))
    sorted_as = sorted(all_AS_values, reverse=True)
    if len(sorted_as) < 2:
        return 60
    gap = sorted_as[0] - sorted_as[1]
    if gap >= 30: return 60
    if gap >= 20: return 50
    if gap >= 10: return 40
    if gap >= 5: return 25
    if gap >= 1: return 10
    return 3


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", help="Bhmem BAM")
    ap.add_argument("bisulfite_genome", nargs="?", default=DEFAULT_BISULFITE)
    ap.add_argument("--bwa", default=DEFAULT_BWA)
    ap.add_argument("--max-reads", type=int, default=1000)
    ap.add_argument("-o", "--output-prefix", default="synthetic_mapq_batch")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    ct_fa_path = os.path.join(args.bisulfite_genome, "CT_conversion/genome_mfa.CT_conversion.fa")
    ga_fa_path = os.path.join(args.bisulfite_genome, "GA_conversion/genome_mfa.GA_conversion.fa")
    for fa in (ct_fa_path, ga_fa_path):
        if not os.path.isfile(fa + ".bwt"):
            print(f"ERROR: BWA index missing: {fa}.bwt", file=sys.stderr)
            sys.exit(1)

    ct_fa = pysam.FastaFile(ct_fa_path)
    ga_fa = pysam.FastaFile(ga_fa_path)
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Step 1: Extract reads and build converted query FASTQs
    print("Step 1: Extracting reads from BAM...", flush=True)
    reads_info = {}  # idx -> {name, is_r2, nm_tag, mapq, chrom, pos, cigar_tuples, seq}
    # 4 trials: (conv_label, query_conversion, index_fa)
    # v1 = direct conversion, v2 = RC(conv(RC(seq)))
    fq_data = {"CT_v1": [], "GA_v1": [], "CT_v2": [], "GA_v2": []}

    n = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None or not read.has_tag("NM"):
            continue
        n += 1
        if n > args.max_reads:
            break

        idx = n - 1
        seq = read.query_sequence.upper()
        is_r2 = bool(read.is_paired and read.is_read2)

        reads_info[idx] = {
            "name": read.query_name,
            "is_r2": is_r2,
            "nm_tag": int(read.get_tag("NM")),
            "mapq": int(read.mapping_quality),
            "chrom": read.reference_name,
            "pos": read.reference_start,
            "cigar_tuples": read.cigartuples,
            "seq": seq,
        }

        v1 = _convert_query_pbat(seq, is_r2)
        v2 = _reverse_complement_dna(_convert_query_pbat(_reverse_complement_dna(seq), is_r2))

        read_id = f"r{idx}"
        qual = "I" * len(v1)
        for label, qseq in [("CT_v1", v1), ("GA_v1", v1), ("CT_v2", v2), ("GA_v2", v2)]:
            fq_data[label].append(f"@{read_id}\n{qseq}\n+\n{qual}\n")

    bam.close()
    total_reads = len(reads_info)
    print(f"  {total_reads} reads extracted", flush=True)

    # Step 2: Write FASTQs and run bwa mem -a
    print("Step 2: Running bwa mem -a (4 batches)...", flush=True)
    trial_hits = {}  # label -> {read_idx -> list of (chrom, pos, AS)}

    tmpdir = tempfile.mkdtemp(prefix="synth_mapq_")
    trial_configs = [
        ("CT_v1", ct_fa_path),
        ("GA_v1", ga_fa_path),
        ("CT_v2", ct_fa_path),
        ("GA_v2", ga_fa_path),
    ]

    for label, idx_fa in trial_configs:
        fq_path = os.path.join(tmpdir, f"{label}.fq")
        sam_path = os.path.join(tmpdir, f"{label}.sam")
        with open(fq_path, "w") as f:
            f.writelines(fq_data[label])

        print(f"  {label}: aligning against {os.path.basename(idx_fa)}...", flush=True)
        with open(sam_path, "w") as out:
            subprocess.run(
                [args.bwa, "mem", "-a", "-t", str(args.threads), idx_fa, fq_path],
                stdout=out, stderr=subprocess.DEVNULL, check=True,
            )

        # Parse hits
        hits = defaultdict(list)
        with open(sam_path) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                fields = line.split("\t")
                if len(fields) < 11:
                    continue
                flag = int(fields[1])
                if flag & 4:
                    continue
                rname = fields[0]
                if not rname.startswith("r"):
                    continue
                ridx = int(rname[1:])
                chrom = fields[2].replace("_CT_converted", "").replace("_GA_converted", "")
                pos = int(fields[3]) - 1
                as_val = None
                for tag in fields[11:]:
                    if tag.startswith("AS:i:"):
                        as_val = int(tag.split(":")[2])
                        break
                if as_val is not None:
                    hits[ridx].append({"chrom": chrom, "pos": pos, "AS": as_val})

        trial_hits[label] = dict(hits)
        os.unlink(fq_path)
        os.unlink(sam_path)

    os.rmdir(tmpdir)
    print("  BWA done", flush=True)

    # Step 3: For each read, compute per-trial synthetic MAPQ and pick
    print("Step 3: Computing synthetic MAPQ and picking trials...", flush=True)
    correct_synth = 0
    correct_min = 0
    results = []

    for idx, info in reads_info.items():
        seq = info["seq"]
        is_r2 = info["is_r2"]
        v1 = _convert_query_pbat(seq, is_r2)
        v2 = _reverse_complement_dna(_convert_query_pbat(_reverse_complement_dna(seq), is_r2))

        trials = []
        for label, conv_q in [("CT_v1", v1), ("GA_v1", v1), ("CT_v2", v2), ("GA_v2", v2)]:
            genome = label[:2]
            fa = ct_fa if genome == "CT" else ga_fa
            c_name = bisulfite_converted_contig_name(fa, info["chrom"], genome)
            if c_name is None:
                continue

            result = _bwa_score_at_pos(conv_q, fa, c_name, info["pos"], info["cigar_tuples"])
            if result is None:
                continue
            as_fixed, nm_fixed = result

            # Get all AS from BWA hits for this trial
            hits = trial_hits.get(label, {}).get(idx, [])
            all_as = [h["AS"] for h in hits]

            synth_mq = _estimate_mapq(as_fixed, all_as)

            trials.append({
                "label": label,
                "AS": as_fixed,
                "NM": nm_fixed,
                "synth_mapq": synth_mq,
                "n_hits": len(hits),
            })

        if not trials:
            continue

        # Pick by Bhmem cascade
        best_synth = max(trials, key=lambda t: (t["synth_mapq"], t["AS"], -t["NM"]))
        best_min = min(trials, key=lambda t: t["NM"])

        nm_tag = info["nm_tag"]
        sc = best_synth["NM"] == nm_tag
        mc = best_min["NM"] == nm_tag
        if sc: correct_synth += 1
        if mc: correct_min += 1

        results.append({
            "name": info["name"],
            "is_r2": info["is_r2"],
            "nm_tag": nm_tag,
            "bam_mapq": info["mapq"],
            "synth_nm": best_synth["NM"],
            "synth_mq": best_synth["synth_mapq"],
            "synth_label": best_synth["label"],
            "min_nm": best_min["NM"],
            "synth_correct": sc,
            "min_correct": mc,
            "all_trials": trials,
        })

    ct_fa.close()
    ga_fa.close()

    # Report
    total = len(results)
    print(f"\n{'='*60}")
    print(f"Synthetic MAPQ trial pick ({total} reads)")
    print(f"{'='*60}\n")

    print(f"Synthetic MAPQ pick matches NM tag: {correct_synth}/{total} ({100*correct_synth/total:.1f}%)")
    print(f"min(trials) matches NM tag:         {correct_min}/{total} ({100*correct_min/total:.1f}%)")

    for mate, flag in [("R1", False), ("R2", True)]:
        sub = [r for r in results if r["is_r2"] == flag]
        if not sub:
            continue
        ns = len(sub)
        cs = sum(1 for r in sub if r["synth_correct"])
        cm = sum(1 for r in sub if r["min_correct"])
        print(f"\n  {mate} ({ns} reads):")
        print(f"    Synthetic MAPQ: {cs}/{ns} ({100*cs/ns:.1f}%)")
        print(f"    min(trials):    {cm}/{ns} ({100*cm/ns:.1f}%)")

        # Where synth != min and synth is correct
        synth_wins = sum(1 for r in sub if r["synth_correct"] and not r["min_correct"])
        min_wins = sum(1 for r in sub if r["min_correct"] and not r["synth_correct"])
        both_wrong = sum(1 for r in sub if not r["synth_correct"] and not r["min_correct"])
        both_right = sum(1 for r in sub if r["synth_correct"] and r["min_correct"])
        print(f"    Both correct:   {both_right}/{ns} ({100*both_right/ns:.1f}%)")
        print(f"    Synth only:     {synth_wins}/{ns} ({100*synth_wins/ns:.1f}%)")
        print(f"    Min only:       {min_wins}/{ns} ({100*min_wins/ns:.1f}%)")
        print(f"    Both wrong:     {both_wrong}/{ns} ({100*both_wrong/ns:.1f}%)")

    # TSV
    tsv_path = f"{args.output_prefix}.per_read.tsv"
    with open(tsv_path, "w") as f:
        f.write("read_name\tis_r2\tnm_tag\tbam_mapq\tsynth_nm\tsynth_mapq\tsynth_label\t"
                "min_nm\tsynth_correct\tmin_correct\ttrials\n")
        for r in results:
            trials_str = ";".join(
                f"{t['label']}:MQ={t['synth_mapq']},AS={t['AS']},NM={t['NM']},hits={t['n_hits']}"
                for t in r["all_trials"]
            )
            f.write(f"{r['name']}\t{int(r['is_r2'])}\t{r['nm_tag']}\t{r['bam_mapq']}\t"
                    f"{r['synth_nm']}\t{r['synth_mq']}\t{r['synth_label']}\t"
                    f"{r['min_nm']}\t{int(r['synth_correct'])}\t{int(r['min_correct'])}\t"
                    f"{trials_str}\n")
    print(f"\nWrote {tsv_path}")


if __name__ == "__main__":
    main()
