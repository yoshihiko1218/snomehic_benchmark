#!/usr/bin/env python3
"""Check if bhmem NM:i == yap XR/XG-recomputed NM for reads with same chrom+pos+CIGAR."""

import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bisulfite_corrected_mismatch import (
    bisulfite_converted_contig_name,
    count_nm_style_edit_distance_converted_explicit,
)

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCseq3"
BS = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10_bismark/Bisulfite_Genome"

ct_fa = pysam.FastaFile(os.path.join(BS, "CT_conversion/genome_mfa.CT_conversion.fa"))
ga_fa = pysam.FastaFile(os.path.join(BS, "GA_conversion/genome_mfa.GA_conversion.fa"))


def normalize_yap(qname):
    parts = qname.split("_")
    if len(parts) < 2:
        return None
    return (parts[0], parts[1] == "2")


def apply_conv(seq, c):
    s = seq.upper()
    return s.replace("C", "T") if c == "CT" else s.replace("G", "A")


def yap_nm_fast(read):
    """Compute NM by fetching entire ref span at once (much faster than per-base)."""
    if read.is_unmapped or read.query_sequence is None:
        return -1
    if not read.has_tag("XR") or not read.has_tag("XG"):
        return -1
    xr = str(read.get_tag("XR")).upper()
    xg = str(read.get_tag("XG")).upper()
    fa = ct_fa if xg == "CT" else ga_fa
    rc = bisulfite_converted_contig_name(fa, read.reference_name, xg)
    if rc is None:
        return -1

    qc = apply_conv(read.query_sequence.upper(), xr)
    ct = read.cigartuples
    if not ct:
        return 0

    # Compute ref span length
    ref_len = sum(l for o, l in ct if o in (0, 2, 3, 7, 8))
    ref_start = read.reference_start
    ref_span = fa.fetch(rc, ref_start, ref_start + ref_len).upper()

    q = 0
    r = 0
    subs = 0
    indel_bases = 0

    for op, length in ct:
        if op in (0, 7, 8):  # M, =, X
            for i in range(length):
                if q < len(qc) and r < len(ref_span):
                    if qc[q] != ref_span[r]:
                        subs += 1
                q += 1
                r += 1
        elif op == 1:  # I
            indel_bases += length
            q += length
        elif op in (2, 3):  # D, N
            indel_bases += length
            r += length
        elif op == 4:  # S
            q += length
        elif op == 5:  # H
            pass

    return subs + indel_bases


# Load bhmem — only store (chrom, pos, cigar, nm) keyed by (name, is_r2)
print("Loading bhmem...", flush=True)
bam_b = pysam.AlignmentFile(os.path.join(BASE, "04.bhmem_bam", "SRR21549292.bhmem.bam"), "rb")
# Use compact storage: key -> (chrom_idx, pos, cigar_hash, nm)
# To save memory, store cigar as hash
bhmem_nm = {}       # (name, is_r2) -> nm
bhmem_loc = {}      # (name, is_r2) -> (chrom, pos, cigar)

n = 0
for r in bam_b:
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    if not r.has_tag("NM"):
        continue
    is_r2 = bool(r.is_paired and r.is_read2)
    key = (r.query_name, is_r2)
    bhmem_nm[key] = int(r.get_tag("NM"))
    bhmem_loc[key] = (r.reference_name, r.reference_start, r.cigarstring)
    n += 1

bam_b.close()
print("Loaded %d bhmem reads" % n, flush=True)

# Scan yap
print("Scanning yap...", flush=True)
bam_y = pysam.AlignmentFile(
    os.path.join(BASE, "alignment", "Group22", "bam", "SRR21549292.3C.sorted.bam"), "rb"
)

exact = {"r1": 0, "r2": 0}
total = {"r1": 0, "r2": 0}
examples = []
n_yap = 0
n_shared = 0

for r in bam_y:
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    if r.query_sequence is None:
        continue

    parsed = normalize_yap(r.query_name)
    if parsed is None:
        continue
    if parsed not in bhmem_nm:
        continue

    n_shared += 1
    bchrom, bpos, bcigar = bhmem_loc[parsed]

    if bchrom != r.reference_name:
        continue
    if bpos != r.reference_start:
        continue
    if bcigar != r.cigarstring:
        continue

    # Same chrom + pos + CIGAR — compute yap NM
    ynm = yap_nm_fast(r)
    if ynm < 0:
        continue

    m = "r2" if parsed[1] else "r1"
    total[m] += 1
    bnm = bhmem_nm[parsed]

    if bnm == ynm:
        exact[m] += 1
    elif len(examples) < 20:
        examples.append(
            (parsed[0], m.upper(), bnm, ynm, bcigar[:60], bchrom, bpos,
             str(r.get_tag("XR")), str(r.get_tag("XG")))
        )

    n_yap += 1
    if n_yap % 5000 == 0:
        print("  %d same-cigar reads found..." % n_yap, flush=True)

bam_y.close()
ct_fa.close()
ga_fa.close()

t = total["r1"] + total["r2"]
e = exact["r1"] + exact["r2"]
print("", flush=True)
print("Shared reads scanned: %d" % n_shared)
print("")
print("=== Same chrom + pos + CIGAR ===")
print("Total: %d" % t)
if t:
    print("NM match: %d/%d (%.2f%%)" % (e, t, 100.0 * e / t))
if total["r1"]:
    print("  R1: %d/%d (%.2f%%)" % (exact["r1"], total["r1"], 100.0 * exact["r1"] / total["r1"]))
if total["r2"]:
    print("  R2: %d/%d (%.2f%%)" % (exact["r2"], total["r2"], 100.0 * exact["r2"] / total["r2"]))

if examples:
    print("")
    print("Examples where NM differs (same pos+cigar):")
    for name, mate, bnm, ynm, cig, chrom, pos, xr, xg in examples:
        print(
            "  %s %s: bhmem=%d yap=%d XR=%s XG=%s cig=%s %s:%d"
            % (name, mate, bnm, ynm, xr, xg, cig, chrom, pos)
        )
