#!/usr/bin/env python3
"""
PDF: one page per read, bhmem + yap Bowtie2 + yap Bowtie1 (same sample).
Bowtie1 uses QNAME like SRR21549292.1_1 / .1_2; Bowtie2 uses .1_1_1 / .1_2_2.
"""

import argparse
import subprocess
import sys


def flag_str(flag):
    f = int(flag)
    parts = []
    if f & 1:
        parts.append("paired")
    if f & 2:
        parts.append("proper_pair")
    if f & 4:
        parts.append("unmapped")
    if f & 8:
        parts.append("mate_unmapped")
    if f & 16:
        parts.append("reverse")
    if f & 32:
        parts.append("mate_reverse")
    if f & 64:
        parts.append("R1")
    if f & 128:
        parts.append("R2")
    if f & 256:
        parts.append("secondary")
    if f & 512:
        parts.append("qcfail")
    if f & 1024:
        parts.append("dup")
    if f & 2048:
        parts.append("supplementary")
    return ",".join(parts) if parts else str(f)


def yap_read_end(qname):
    """R1 or R2 from yap QNAME (Bowtie1 or Bowtie2 naming)."""
    parts = qname.split("_")
    if len(parts) < 2:
        return None
    bit = parts[1]
    if not bit:
        return None
    if bit[0] == "1":
        return "R1"
    if bit[0] == "2":
        return "R2"
    return None


def get_r1_r2_info(rows, source):
    r1_chr, r1_pos, r2_chr, r2_pos = None, None, None, None
    if source == "bhmem":
        for line in rows:
            p = line.split("\t")
            if len(p) < 8:
                continue
            flag, rname, pos = int(p[1]), p[2], int(p[3])
            if flag & 64:
                r1_chr, r1_pos = rname, pos
            elif flag & 128:
                r2_chr, r2_pos = rname, pos
    else:
        candidates = {"R1": [], "R2": []}
        for line in rows:
            p = line.split("\t")
            if len(p) < 6:
                continue
            qname, rname, pos, mapq = p[0], p[2], int(p[3]), int(p[4])
            end = yap_read_end(qname)
            if end is None:
                continue
            parts = qname.split("_")
            is_prim = "-" not in parts[-1]
            candidates[end].append((rname, pos, mapq, is_prim))
        for end in ("R1", "R2"):
            if not candidates[end]:
                continue
            prim = [c for c in candidates[end] if c[3]]
            pool = prim if prim else candidates[end]
            best = max(pool, key=lambda c: c[2])
            if end == "R1":
                r1_chr, r1_pos = best[0], best[1]
            else:
                r2_chr, r2_pos = best[0], best[1]
    return r1_chr, r1_pos, r2_chr, r2_pos


def fmt_row(line):
    if not line:
        return ""
    p = line.split("\t")
    if len(p) < 11:
        return line + "\n"
    qname, flag, rname, pos, mapq, cigar = p[0], p[1], p[2], p[3], p[4], p[5]
    rnext, pnext = p[6], p[7] if len(p) > 7 else ""
    seq = p[9][:55] + "..." if len(p[9]) > 55 else p[9]
    nm = ""
    for t in p[11:]:
        if t.startswith("NM:"):
            nm = t.split(":")[-1]
            break
    return (
        f"  QNAME: {qname}\n"
        f"  FLAG: {flag} ({flag_str(flag)})\n"
        f"  RNAME: {rname}  POS: {pos}  MAPQ: {mapq}  CIGAR: {cigar}\n"
        f"  RNEXT: {rnext}  PNEXT: {pnext}  NM: {nm}\n"
        f"  SEQ: {seq}\n"
    )


def dist_str(r1c, r1p, r2c, r2p):
    if not all([r1c, r1p, r2c, r2p]):
        return "N/A"
    if r1c != r2c:
        return "trans (diff chr)"
    return f"{abs(r1p - r2p):,} bp"


def prefetch_yap(proc_stdout, base_ids):
    by_base = {bid: [] for bid in base_ids}
    for line in proc_stdout:
        qname = line.split("\t")[0]
        for bid in base_ids:
            if qname.startswith(bid + "_"):
                by_base[bid].append(line.rstrip())
                break
    return by_base


def main():
    ap = argparse.ArgumentParser(
        description="PDF: bhmem + yap Bowtie2 + yap Bowtie1 for discrepant reads"
    )
    ap.add_argument("subset_tsv", help="yap_high_bhmem_low.tsv")
    ap.add_argument("bhmem_bam", help="bhmem BAM")
    ap.add_argument("yap_bowtie2_bam", help="yap 3C BAM (Bismark Bowtie2)")
    ap.add_argument("yap_bowtie1_bam", help="yap 3C BAM (Bismark Bowtie1)")
    ap.add_argument("-o", "--output", default="discrepant_reads_three.pdf")
    ap.add_argument("-n", "--num", type=int, default=10)
    args = ap.parse_args()

    base_ids = []
    with open(args.subset_tsv) as f:
        next(f)
        for line in f:
            bid = line.split("\t")[0]
            if bid not in base_ids:
                base_ids.append(bid)
            if len(base_ids) >= args.num:
                break

    base_set = set(base_ids)

    proc = subprocess.Popen(
        ["samtools", "view", args.yap_bowtie2_bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    yap2 = prefetch_yap(proc.stdout, base_ids)
    proc.wait()

    proc = subprocess.Popen(
        ["samtools", "view", args.yap_bowtie1_bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    yap1 = prefetch_yap(proc.stdout, base_ids)
    proc.wait()

    bhmem_by_base = {bid: [] for bid in base_ids}
    proc = subprocess.Popen(
        ["samtools", "view", args.bhmem_bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    for line in proc.stdout:
        qname = line.split("\t")[0]
        if qname in base_set:
            bhmem_by_base[qname].append(line.rstrip())
    proc.wait()

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        print("ERROR: matplotlib required", file=sys.stderr)
        sys.exit(1)

    with PdfPages(args.output) as pdf:
        for i, base_id in enumerate(base_ids):
            bh = bhmem_by_base.get(base_id, [])
            y2 = yap2.get(base_id, [])
            y1 = yap1.get(base_id, [])

            b_r = get_r1_r2_info(bh, "bhmem")
            y2_r = get_r1_r2_info(y2, "yap")
            y1_r = get_r1_r2_info(y1, "yap")

            blocks = []
            blocks.append(f"Read {i+1}/{args.num}: {base_id}\n")
            blocks.append("Pipelines: bhmem (BwaMem) | yap Bowtie2 | yap Bowtie1\n")
            blocks.append("--- Summary: R1-R2 distance ---\n")
            blocks.append(f"  bhmem:        {dist_str(*b_r)}\n")
            blocks.append(f"  yap Bowtie2:  {dist_str(*y2_r)}\n")
            blocks.append(f"  yap Bowtie1:  {dist_str(*y1_r)}\n")
            blocks.append(
                "Secondary: check FLAG&256 in raw BAM (not summarized here).\n"
            )
            blocks.append("=== bhmem ===\n")
            blocks.append(f"Rows: {len(bh)}\n")
            for j, row in enumerate(bh):
                blocks.append(f"--- {j+1} ---\n")
                blocks.append(fmt_row(row))
            blocks.append("\n=== yap Bowtie2 (Bismark --bowtie2) ===\n")
            blocks.append(f"Rows: {len(y2)} (QNAME e.g. {base_id}_1_1)\n")
            for j, row in enumerate(y2):
                blocks.append(f"--- {j+1} ---\n")
                blocks.append(fmt_row(row))
            blocks.append("\n=== yap Bowtie1 (Bismark --bowtie1) ===\n")
            blocks.append(f"Rows: {len(y1)} (QNAME e.g. {base_id}_1)\n")
            for j, row in enumerate(y1):
                blocks.append(f"--- {j+1} ---\n")
                blocks.append(fmt_row(row))
            blocks.append("\nSplit suffixes: -l -r -m (m3c-split-reads)\n")

            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis("off")
            ax.text(
                0.02,
                0.98,
                "".join(blocks),
                transform=ax.transAxes,
                fontsize=6.5,
                verticalalignment="top",
                fontfamily="monospace",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.25),
            )
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

    print(f"Wrote {args.output} ({len(base_ids)} pages)")


if __name__ == "__main__":
    main()
