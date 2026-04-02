#!/usr/bin/env python3
"""
Compute % chimeric reads (R1 and R2 on different chromosomes) for the discrepant subset.
Uses RNEXT from bhmem BAM; for yap, looks up both R1 and R2 chr from 3C BAM.
"""

import argparse
import subprocess
import sys
from collections import defaultdict

def _is_r1(flag: int) -> bool:
    return (flag & 64) != 0

def parse_yap_qname(qname: str):
    parts = qname.split("_")
    if len(parts) < 2 or parts[1] not in ("1", "2"):
        return None, None
    return parts[0], (parts[1] == "1")

def load_subset_ids(path: str):
    ids = set()
    with open(path) as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 2:
                ids.add((p[0], bool(int(p[1]))))
    return ids

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("subset_tsv", help="yap_high_bhmem_low.tsv")
    ap.add_argument("bhmem_bam", help="bhmem BAM")
    ap.add_argument("yap_bam", help="yap 3C BAM")
    args = ap.parse_args()

    subset_ids = load_subset_ids(args.subset_tsv)
    subset_base_ids = {k[0] for k in subset_ids}
    print(f"Loaded {len(subset_ids)} reads from {len(subset_base_ids)} base_ids", file=sys.stderr)

    # bhmem: RNEXT is mate's chr; chimeric = (RNAME != RNEXT && RNEXT != '*')
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "4", args.bhmem_bam],
        stdout=subprocess.PIPE, text=True
    )
    bhmem_chimeric = 0
    bhmem_total = 0
    for line in proc.stdout:
        p = line.strip().split("\t")
        if len(p) < 8:
            continue
        qname, flag, rname, rnext = p[0], int(p[1]), p[2], p[6]
        key = (qname, _is_r1(flag))
        if key not in subset_ids:
            continue
        bhmem_total += 1
        if rnext != "*" and rname != rnext:
            bhmem_chimeric += 1
    proc.wait()

    # yap: need both R1 and R2 chr per base_id; stream and collect for subset base_ids
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "4", args.yap_bam],
        stdout=subprocess.PIPE, text=True
    )
    chr_by_base = defaultdict(dict)  # base_id -> {is_r1: (chr, mapq)}
    for line in proc.stdout:
        p = line.strip().split("\t")
        if len(p) < 5:
            continue
        base, is_r1 = parse_yap_qname(p[0])
        if base is None or base not in subset_base_ids:
            continue
        rname, mapq = p[2], int(p[4])
        if is_r1 not in chr_by_base[base] or mapq > chr_by_base[base][is_r1][1]:
            chr_by_base[base][is_r1] = (rname, mapq)
    proc.wait()

    # For reads in subset, get chimeric from chr_r1 != chr_r2
    yap_chimeric = 0
    yap_total = 0
    for (base_id, is_r1) in subset_ids:
        if base_id not in chr_by_base or True not in chr_by_base[base_id] or False not in chr_by_base[base_id]:
            continue
        c1 = chr_by_base[base_id][True][0]
        c2 = chr_by_base[base_id][False][0]
        yap_total += 1
        if c1 != c2:
            yap_chimeric += 1

    print("\n--- Chimeric reads (R1 and R2 on different chromosomes) ---")
    print(f"bhmem: {bhmem_chimeric} / {bhmem_total} = {100*bhmem_chimeric/bhmem_total:.1f}%" if bhmem_total else "bhmem: N/A")
    print(f"yap:   {yap_chimeric} / {yap_total} = {100*yap_chimeric/yap_total:.1f}%" if yap_total else "yap: N/A")

if __name__ == "__main__":
    main()
