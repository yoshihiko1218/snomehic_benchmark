#!/usr/bin/env python3
"""
Create PDF with one discrepant read per page (bhmem low MAPQ, yap high MAPQ).
Shows ALL rows from both BAMs including any secondary/split alignments.
"""

import argparse
import subprocess
import sys
from io import StringIO

def flag_str(flag):
    f = int(flag)
    parts = []
    if f & 1: parts.append("paired")
    if f & 2: parts.append("proper_pair")
    if f & 4: parts.append("unmapped")
    if f & 8: parts.append("mate_unmapped")
    if f & 16: parts.append("reverse")
    if f & 32: parts.append("mate_reverse")
    if f & 64: parts.append("R1")
    if f & 128: parts.append("R2")
    if f & 256: parts.append("secondary")
    if f & 512: parts.append("qcfail")
    if f & 1024: parts.append("dup")
    if f & 2048: parts.append("supplementary")
    return ",".join(parts) if parts else str(f)

def get_r1_r2_info(rows, source):
    """Extract R1 and R2 (chr, pos) for distance calc. source='bhmem' or 'yap'."""
    r1_chr, r1_pos, r2_chr, r2_pos = None, None, None, None
    if source == "bhmem":
        for line in rows:
            p = line.split("\t")
            if len(p) < 8:
                continue
            flag, rname, pos = int(p[1]), p[2], int(p[3])
            if flag & 64:  # R1
                r1_chr, r1_pos = rname, pos
            elif flag & 128:  # R2
                r2_chr, r2_pos = rname, pos
    else:  # yap: _1_ = R1, _2_ = R2; prefer primary (no -l/-r/-m)
        candidates = {"R1": [], "R2": []}
        for line in rows:
            p = line.split("\t")
            if len(p) < 6:
                continue
            qname, rname, pos, mapq = p[0], p[2], int(p[3]), int(p[4])
            parts = qname.split("_")
            if len(parts) < 2:
                continue
            end = "R1" if parts[1] == "1" else "R2"
            is_prim = "-" not in parts[-1]
            candidates[end].append((rname, pos, mapq, is_prim))
        for end in ("R1", "R2"):
            if not candidates[end]:
                continue
            prim = [c for c in candidates[end] if c[3]]
            pool = prim if prim else candidates[end]
            best = max(pool, key=lambda c: c[2])  # max MAPQ
            if end == "R1":
                r1_chr, r1_pos = best[0], best[1]
            else:
                r2_chr, r2_pos = best[0], best[1]
    return r1_chr, r1_pos, r2_chr, r2_pos

def fmt_row(line, source):
    if not line:
        return ""
    p = line.split("\t")
    if len(p) < 11:
        return line + "\n"
    qname, flag, rname, pos, mapq, cigar = p[0], p[1], p[2], p[3], p[4], p[5]
    rnext, pnext = p[6], p[7] if len(p) > 7 else ""
    seq = p[9][:60] + "..." if len(p[9]) > 60 else p[9]
    tags = "\t".join(p[11:])
    nm = ""
    for t in tags.split("\t"):
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("subset_tsv", help="yap_high_bhmem_low.tsv")
    ap.add_argument("bhmem_bam", help="bhmem BAM")
    ap.add_argument("yap_bam", help="yap 3C BAM")
    ap.add_argument("-o", "--output", default="discrepant_reads.pdf", help="Output PDF")
    ap.add_argument("-n", "--num", type=int, default=10, help="Number of reads (base_ids)")
    args = ap.parse_args()

    # Load subset and pick base_ids (unique)
    base_ids = []
    with open(args.subset_tsv) as f:
        next(f)
        for line in f:
            bid = line.split("\t")[0]
            if bid not in base_ids:
                base_ids.append(bid)
            if len(base_ids) >= args.num:
                break

    # Pre-fetch yap rows for all base_ids in one pass
    yap_by_base = {bid: [] for bid in base_ids}
    proc = subprocess.Popen(
        ["samtools", "view", args.yap_bam],
        stdout=subprocess.PIPE, text=True
    )
    for line in proc.stdout:
        qname = line.split("\t")[0]
        for bid in base_ids:
            if qname.startswith(bid + "_"):
                yap_by_base[bid].append(line.rstrip())
                break
    proc.wait()

    # Pre-fetch bhmem rows for all base_ids in one pass
    bhmem_by_base = {bid: [] for bid in base_ids}
    proc = subprocess.Popen(
        ["samtools", "view", args.bhmem_bam],
        stdout=subprocess.PIPE, text=True
    )
    base_set = set(base_ids)
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
        print("ERROR: matplotlib required for PDF", file=sys.stderr)
        sys.exit(1)

    with PdfPages(args.output) as pdf:
        for i, base_id in enumerate(base_ids):
            bhmem_rows = bhmem_by_base.get(base_id, [])
            yap_rows = yap_by_base.get(base_id, [])

            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis("off")

            # Summary: R1-R2 distance
            b_r1c, b_r1p, b_r2c, b_r2p = get_r1_r2_info(bhmem_rows, "bhmem")
            y_r1c, y_r1p, y_r2c, y_r2p = get_r1_r2_info(yap_rows, "yap")

            def dist_str(r1c, r1p, r2c, r2p):
                if not all([r1c, r1p, r2c, r2p]):
                    return "N/A"
                if r1c != r2c:
                    return "trans (diff chr)"
                return f"{abs(r1p - r2p):,} bp"

            b_dist = dist_str(b_r1c, b_r1p, b_r2c, b_r2p)
            y_dist = dist_str(y_r1c, y_r1p, y_r2c, y_r2p)

            blocks = []
            blocks.append(f"Read {i+1}/{args.num}: {base_id}\n")
            blocks.append("--- Summary: R1-R2 distance ---\n")
            blocks.append(f"  bhmem: {b_dist}\n")
            blocks.append(f"  yap:   {y_dist}\n")
            blocks.append("Secondary alignments: Neither BAM has FLAG&256 (secondary).\n")
            blocks.append("=== bhmem BAM ===\n")
            blocks.append(f"Rows: {len(bhmem_rows)}\n")
            for j, row in enumerate(bhmem_rows):
                blocks.append(f"--- Row {j+1} ---\n")
                blocks.append(fmt_row(row, "bhmem"))
            blocks.append("\n=== yap 3C BAM ===\n")
            blocks.append(f"Rows: {len(yap_rows)} (full + split with -l/-r/-m)\n")
            for j, row in enumerate(yap_rows):
                blocks.append(f"--- Row {j+1} ---\n")
                blocks.append(fmt_row(row, "yap"))
            blocks.append("\n(Split reads: -l left, -r right, -m middle from m3c-split-reads)")

            full_text = "".join(blocks)
            ax.text(0.02, 0.98, full_text, transform=ax.transAxes, fontsize=8,
                    verticalalignment="top", fontfamily="monospace",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.3))

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

    print(f"Wrote {args.output} ({len(base_ids)} pages)")

if __name__ == "__main__":
    main()
