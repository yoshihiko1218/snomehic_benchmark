#!/usr/bin/env python3
"""
Extract read ID and MAPQ from BAM files for MAPQ comparison between yap (Bowtie2) and bhmem (BwaMem).

Read ID mapping:
- bhmem: QNAME = SRR21549292.N (same for R1/R2 in a pair; use FLAG to distinguish)
- yap:   QNAME = SRR21549292.N_1_1 or SRR21549292.N_2_2 (N_1 = R1, N_2 = R2)

Usage:
  python extract_mapq.py bhmem /path/to/sample.bhmem.bam -o bhmem_mapq.tsv
  python extract_mapq.py yap   /path/to/sample.3C.sorted.bam -o yap_mapq.tsv
"""

import argparse
import subprocess
import sys
from collections import defaultdict

USE_PYSAM = False
try:
    import pysam
    USE_PYSAM = True
except ImportError:
    pass


def is_r1(flag: int) -> bool:
    """True if read is first in pair (R1)."""
    return (flag & 64) != 0


def parse_yap_qname(qname: str):
    """
    Parse yap QNAME to (base_id, is_r1).
    e.g. SRR21549292.1_1_1 -> (SRR21549292.1, True)
         SRR21549292.1_2_2 -> (SRR21549292.1, False)
         SRR21549292.3_2_2-l -> (SRR21549292.3, False)  [split read]
    """
    parts = qname.split("_")
    if len(parts) < 2:
        return None, None
    base = parts[0]  # SRR21549292.1
    strand = parts[1]  # "1" or "2"
    if strand not in ("1", "2"):
        return None, None
    is_r1 = strand == "1"
    return base, is_r1


def extract_bhmem_samtools(bam_path: str, out_path: str, min_mapq: int = -1) -> int:
    """Extract using samtools view (no pysam)."""
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "4", bam_path],
        stdout=subprocess.PIPE, text=True
    )
    count = 0
    with open(out_path, "w") as out:
        out.write("base_id\tis_r1\tmapq\n")
        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t", 11)
            if len(parts) < 5:
                continue
            qname, flag, _, _, mapq = parts[0], int(parts[1]), parts[2], parts[3], int(parts[4])
            if min_mapq >= 0 and mapq < min_mapq:
                continue
            out.write(f"{qname}\t{int(is_r1(flag))}\t{mapq}\n")
            count += 1
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools view failed with code {proc.returncode}")
    return count


def extract_bhmem_pysam(bam_path: str, out_path: str, min_mapq: int = -1) -> int:
    """Extract using pysam."""
    count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam, open(out_path, "w") as out:
        out.write("base_id\tis_r1\tmapq\n")
        for r in bam:
            if r.is_unmapped:
                continue
            if min_mapq >= 0 and r.mapping_quality < min_mapq:
                continue
            base_id = r.query_name
            out.write(f"{base_id}\t{int(is_r1(r.flag))}\t{r.mapping_quality}\n")
            count += 1
    return count


def extract_bhmem(bam_path: str, out_path: str, min_mapq: int = -1) -> int:
    if USE_PYSAM:
        return extract_bhmem_pysam(bam_path, out_path, min_mapq)
    return extract_bhmem_samtools(bam_path, out_path, min_mapq)


def extract_yap_samtools(bam_path: str, out_path: str, min_mapq: int = -1, primary_only: bool = False) -> int:
    """Extract using samtools view (no pysam)."""
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "4", bam_path],
        stdout=subprocess.PIPE, text=True
    )
    by_key = defaultdict(list)
    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t", 11)
        if len(parts) < 5:
            continue
        qname, _, _, _, mapq = parts[0], parts[1], parts[2], parts[3], int(parts[4])
        if min_mapq >= 0 and mapq < min_mapq:
            continue
        base, is_r1 = parse_yap_qname(qname)
        if base is None:
            continue
        last_part = qname.split("_")[-1] if "_" in qname else ""
        is_primary = "-" not in last_part
        by_key[(base, is_r1)].append((mapq, is_primary))
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools view failed with code {proc.returncode}")

    with open(out_path, "w") as out:
        out.write("base_id\tis_r1\tmapq\n")
        count = 0
        for (base_id, is_r1), candidates in sorted(by_key.items()):
            if primary_only:
                primaries = [mq for mq, prim in candidates if prim]
                mq = max(primaries) if primaries else max(mq for mq, _ in candidates)
            else:
                mq = max(mq for mq, _ in candidates)
            out.write(f"{base_id}\t{int(is_r1)}\t{mq}\n")
            count += 1
    return count


def extract_yap_pysam(bam_path: str, out_path: str, min_mapq: int = -1, primary_only: bool = False) -> int:
    """Extract using pysam."""
    by_key = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped:
                continue
            if min_mapq >= 0 and r.mapping_quality < min_mapq:
                continue
            base, is_r1 = parse_yap_qname(r.query_name)
            if base is None:
                continue
            last_part = r.query_name.split("_")[-1] if "_" in r.query_name else ""
            is_primary = "-" not in last_part
            by_key[(base, is_r1)].append((r.mapping_quality, is_primary))

    with open(out_path, "w") as out:
        out.write("base_id\tis_r1\tmapq\n")
        count = 0
        for (base_id, is_r1), candidates in sorted(by_key.items()):
            if primary_only:
                primaries = [mq for mq, prim in candidates if prim]
                mq = max(primaries) if primaries else max(mq for mq, _ in candidates)
            else:
                mq = max(mq for mq, _ in candidates)
            out.write(f"{base_id}\t{int(is_r1)}\t{mq}\n")
            count += 1
    return count


def extract_yap(bam_path: str, out_path: str, min_mapq: int = -1, primary_only: bool = False) -> int:
    if USE_PYSAM:
        return extract_yap_pysam(bam_path, out_path, min_mapq, primary_only)
    return extract_yap_samtools(bam_path, out_path, min_mapq, primary_only)


def main():
    ap = argparse.ArgumentParser(description="Extract MAPQ from bhmem or yap BAM")
    ap.add_argument("mode", choices=("bhmem", "yap"), help="BAM source: bhmem or yap")
    ap.add_argument("bam", help="Path to BAM file")
    ap.add_argument("-o", "--output", required=True, help="Output TSV path")
    ap.add_argument("--min-mapq", type=int, default=-1, help="Skip reads with MAPQ < this")
    ap.add_argument("--primary-only", action="store_true", help="[yap only] Use primary alignments only (exclude split reads)")
    args = ap.parse_args()

    if args.mode == "bhmem":
        n = extract_bhmem(args.bam, args.output, args.min_mapq)
    else:
        n = extract_yap(args.bam, args.output, args.min_mapq, args.primary_only)

    print(f"Wrote {n} records to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
