#!/usr/bin/env python3
"""
For reads where yap MAPQ > 30 and bhmem MAPQ < 30, extract alignment details
from both BAMs and summarize: same location, distance, mismatches, etc.
"""

import argparse
import re
import subprocess
import sys
from collections import defaultdict

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def is_r1(flag: int) -> bool:
    return (flag & 64) != 0


def parse_yap_qname(qname: str):
    parts = qname.split("_")
    if len(parts) < 2:
        return None, None
    base, strand = parts[0], parts[1]
    if strand not in ("1", "2"):
        return None, None
    return base, (strand == "1")


def parse_tag(tags_str: str, key: str):
    """Parse SAM tags like 'NM:i:5' -> 5."""
    for t in tags_str.split("\t"):
        if t.startswith(key + ":"):
            if "i:" in t:
                return int(t.split(":")[-1])
            return t.split(":")[-1]
    return None


def load_subset_ids(path: str):
    """Load (base_id, is_r1) from yap_high_bhmem_low.tsv."""
    ids = set()
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                base_id, is_r1 = parts[0], int(parts[1])
                ids.add((base_id, bool(is_r1)))
    return ids


def extract_bam_info_yap(bam_path: str, subset_ids: set):
    """Yap: may have multiple alignments per (base_id, is_r1); take primary or best MAPQ."""
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "4", bam_path],
        stdout=subprocess.PIPE, text=True
    )
    by_key = defaultdict(list)
    for line in proc.stdout:
        parts = line.strip().split("\t", 11)
        if len(parts) < 11:
            continue
        qname, flag, rname, pos, mapq = parts[0], int(parts[1]), parts[2], int(parts[3]), int(parts[4])
        seq = parts[9]
        tags = parts[11] if len(parts) > 11 else ""

        base, is_r1 = parse_yap_qname(qname)
        if base is None:
            continue
        key = (base, is_r1)
        if key not in subset_ids:
            continue

        nm = parse_tag(tags, "NM")
        strand = "-" if (flag & 16) else "+"
        is_primary = "-" not in (qname.split("_")[-1] if "_" in qname else "")

        by_key[key].append({
            "chr": rname, "pos": pos, "nm": nm if nm is not None else -1,
            "seq_len": len(seq), "strand": strand, "mapq": mapq, "primary": is_primary,
        })
    proc.wait()

    info = {}
    for key, candidates in by_key.items():
        primaries = [c for c in candidates if c["primary"]]
        pool = primaries if primaries else candidates
        best = max(pool, key=lambda c: c["mapq"])
        info[key] = {k: v for k, v in best.items() if k != "primary"}
    return info


def extract_bam_info_bhmem(bam_path: str, subset_ids: set):
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "4", bam_path],
        stdout=subprocess.PIPE, text=True
    )
    info = {}
    for line in proc.stdout:
        parts = line.strip().split("\t", 11)
        if len(parts) < 11:
            continue
        qname, flag, rname, pos, mapq = parts[0], int(parts[1]), parts[2], int(parts[3]), int(parts[4])
        seq = parts[9]
        tags = parts[11] if len(parts) > 11 else ""

        base_id = qname
        is_r1_flag = is_r1(flag)
        key = (base_id, is_r1_flag)
        if key not in subset_ids:
            continue

        nm = parse_tag(tags, "NM")
        strand = "-" if (flag & 16) else "+"
        info[key] = {
            "chr": rname, "pos": pos, "nm": nm if nm is not None else -1,
            "seq_len": len(seq), "strand": strand, "mapq": mapq,
        }
    proc.wait()
    return info


def main():
    ap = argparse.ArgumentParser(description="Summarize reads where yap MAPQ>30 and bhmem MAPQ<30")
    ap.add_argument("subset_tsv", help="yap_high_bhmem_low.tsv (base_id, is_r1, mapq_bhmem, mapq_yap)")
    ap.add_argument("bhmem_bam", help="bhmem BAM")
    ap.add_argument("yap_bam", help="yap 3C BAM")
    ap.add_argument("-o", "--output-prefix", default="discrepant_summary", help="Output prefix")
    ap.add_argument("--mapq-tsv", help="Optional: yap_high_bhmem_low.tsv to add mapq_bhmem, mapq_yap to detail")
    args = ap.parse_args()

    subset_ids = load_subset_ids(args.subset_tsv)
    print(f"Loaded {len(subset_ids)} read keys from subset", file=sys.stderr)

    bhmem_info = extract_bam_info_bhmem(args.bhmem_bam, subset_ids)
    print(f"Found {len(bhmem_info)} in bhmem BAM", file=sys.stderr)

    yap_info = extract_bam_info_yap(args.yap_bam, subset_ids)
    print(f"Found {len(yap_info)} in yap BAM", file=sys.stderr)

    keys = sorted(set(bhmem_info) & set(yap_info))
    print(f"Overlap: {len(keys)} reads with info in both", file=sys.stderr)

    # Compute per-read stats
    same_chr = []
    same_pos = []
    pos_dist = []
    nm_bhmem = []
    nm_yap = []
    nm_diff = []
    same_strand = []
    seq_lens = []

    rows = []
    for k in keys:
        b, y = bhmem_info[k], yap_info[k]
        sc = b["chr"] == y["chr"]
        same_chr.append(sc)
        if sc:
            d = abs(b["pos"] - y["pos"])
            pos_dist.append(d)
            same_pos.append(d == 0)
        else:
            pos_dist.append(None)
            same_pos.append(False)
        nm_bhmem.append(b["nm"] if b["nm"] >= 0 else None)
        nm_yap.append(y["nm"] if y["nm"] >= 0 else None)
        if b["nm"] >= 0 and y["nm"] >= 0:
            nm_diff.append(b["nm"] - y["nm"])
        same_strand.append(b["strand"] == y["strand"])
        seq_lens.append(b["seq_len"])

        rows.append({
            "base_id": k[0], "is_r1": k[1],
            "chr_bhmem": b["chr"], "pos_bhmem": b["pos"], "nm_bhmem": b["nm"],
            "chr_yap": y["chr"], "pos_yap": y["pos"], "nm_yap": y["nm"],
            "same_chr": sc, "pos_dist": abs(b["pos"] - y["pos"]) if sc else None,
            "same_strand": b["strand"] == y["strand"],
        })

    # Summary stats
    n = len(keys)
    n_same_chr = sum(same_chr)
    n_same_pos = sum(same_pos)
    valid_dist = [d for d in pos_dist if d is not None]
    valid_nm_b = [x for x in nm_bhmem if x is not None]
    valid_nm_y = [x for x in nm_yap if x is not None]

    stats = []
    stats.append(("N_reads", n))
    stats.append(("N_same_chr", n_same_chr))
    stats.append(("N_diff_chr", n - n_same_chr))
    stats.append(("pct_same_chr", 100 * n_same_chr / n if n else 0))
    stats.append(("N_same_pos", n_same_pos))
    stats.append(("pct_same_pos", 100 * n_same_pos / n if n else 0))
    stats.append(("N_same_strand", sum(same_strand)))
    stats.append(("pct_same_strand", 100 * sum(same_strand) / n if n else 0))

    if valid_dist:
        if HAS_NUMPY:
            stats.append(("pos_dist_mean", float(np.mean(valid_dist))))
            stats.append(("pos_dist_median", float(np.median(valid_dist))))
            stats.append(("pos_dist_max", int(np.max(valid_dist))))
            stats.append(("pos_dist_min", int(np.min(valid_dist))))
        else:
            stats.append(("pos_dist_mean", sum(valid_dist) / len(valid_dist)))
            stats.append(("pos_dist_median", sorted(valid_dist)[len(valid_dist) // 2]))
            stats.append(("pos_dist_max", max(valid_dist)))
            stats.append(("pos_dist_min", min(valid_dist)))

        # Distance bins
        d0 = sum(1 for d in valid_dist if d == 0)
        d1_100 = sum(1 for d in valid_dist if 1 <= d <= 100)
        d101_1k = sum(1 for d in valid_dist if 101 <= d <= 1000)
        d1k_10k = sum(1 for d in valid_dist if 1001 <= d <= 10000)
        d10k_plus = sum(1 for d in valid_dist if d > 10000)
        stats.append(("pos_dist_eq_0", d0))
        stats.append(("pos_dist_1_100bp", d1_100))
        stats.append(("pos_dist_101_1kb", d101_1k))
        stats.append(("pos_dist_1kb_10kb", d1k_10k))
        stats.append(("pos_dist_gt_10kb", d10k_plus))

    if valid_nm_b:
        if HAS_NUMPY:
            stats.append(("nm_bhmem_mean", float(np.mean(valid_nm_b))))
            stats.append(("nm_bhmem_median", float(np.median(valid_nm_b))))
            stats.append(("nm_yap_mean", float(np.mean(valid_nm_y))))
            stats.append(("nm_yap_median", float(np.median(valid_nm_y))))
        else:
            stats.append(("nm_bhmem_mean", sum(valid_nm_b) / len(valid_nm_b)))
            stats.append(("nm_bhmem_median", sorted(valid_nm_b)[len(valid_nm_b) // 2]))
            stats.append(("nm_yap_mean", sum(valid_nm_y) / len(valid_nm_y)))
            stats.append(("nm_yap_median", sorted(valid_nm_y)[len(valid_nm_y) // 2]))

    if nm_diff and HAS_NUMPY:
        stats.append(("nm_diff_mean", float(np.mean(nm_diff))))
        stats.append(("nm_diff_median", float(np.median(nm_diff))))

    if seq_lens and HAS_NUMPY:
        stats.append(("seq_len_mean", float(np.mean(seq_lens))))
        stats.append(("seq_len_median", float(np.median(seq_lens))))

    # Write summary
    summary_path = f"{args.output_prefix}.summary.txt"
    with open(summary_path, "w") as f:
        f.write("# Reads where yap MAPQ>30 and bhmem MAPQ<30\n")
        f.write("# NM note: bhmem aligns to bisulfite-converted ref (low NM); yap/Bismark uses unconverted ref (C-to-T count as mismatches, higher NM)\n")
        for name, val in stats:
            f.write(f"{name}\t{val}\n")
    print(f"Wrote summary to {summary_path}")

    # Optional: load MAPQ from subset TSV
    mapq_by_key = {}
    if args.mapq_tsv:
        with open(args.mapq_tsv) as f:
            next(f)
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 4:
                    mapq_by_key[(p[0], bool(int(p[1])))] = (int(p[2]), int(p[3]))

    # Write per-read detail
    detail_path = f"{args.output_prefix}.detail.tsv"
    with open(detail_path, "w") as f:
        hdr = "base_id\tis_r1\tchr_bhmem\tpos_bhmem\tchr_yap\tpos_yap\tsame_chr\tpos_dist\tsame_strand\tnm_bhmem\tnm_yap"
        if mapq_by_key:
            hdr += "\tmapq_bhmem\tmapq_yap"
        f.write(hdr + "\n")
        for r in rows:
            k = (r["base_id"], r["is_r1"])
            pd = "" if r["pos_dist"] is None else r["pos_dist"]
            line = f"{r['base_id']}\t{r['is_r1']}\t{r['chr_bhmem']}\t{r['pos_bhmem']}\t{r['chr_yap']}\t{r['pos_yap']}\t{r['same_chr']}\t{pd}\t{r['same_strand']}\t{r['nm_bhmem']}\t{r['nm_yap']}"
            if k in mapq_by_key:
                mb, my = mapq_by_key[k]
                line += f"\t{mb}\t{my}"
            f.write(line + "\n")
    print(f"Wrote detail to {detail_path}")

    # Print summary to stdout
    print("\n--- Summary ---")
    for name, val in stats:
        print(f"  {name}: {val}")


if __name__ == "__main__":
    main()
