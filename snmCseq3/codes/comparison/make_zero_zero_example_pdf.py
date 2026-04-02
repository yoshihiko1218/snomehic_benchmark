#!/usr/bin/env python3
"""
Example PDF for zero-zero reads (NM_bhmem=0, yap corrected mismatch=0): pick K reads
per alignment-match category from zero_zero_location_compare.tsv, one page per read
end, same layout as make_discrepant_pdf.py (bhmem + yap 3C BAM, full fragment rows).

Requires discrepant_mismatch_report*.per_read.tsv to filter yap MAPQ > 30 (original
cohort) and to print NM / bisulfite-corrected mismatch stats. Optionally requires both
R1 and R2 alignments present so each page shows a full paired fragment.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("ERROR: matplotlib required", file=sys.stderr)
    sys.exit(1)


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


def parse_yap_read_end(qname: str) -> bool | None:
    """True=R1, False=R2, None=unknown."""
    parts = qname.split("_")
    if len(parts) < 2:
        return None
    bit = parts[1]
    if bit.startswith("1"):
        return True
    if bit.startswith("2"):
        return False
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
            end = parse_yap_read_end(qname)
            if end is None:
                continue
            key = "R1" if end else "R2"
            parts = qname.split("_")
            is_prim = "-" not in parts[-1]
            candidates[key].append((rname, pos, mapq, is_prim))
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


def dist_str(r1c, r1p, r2c, r2p):
    if not all([r1c, r1p, r2c, r2p]):
        return "N/A"
    if r1c != r2c:
        return "trans (diff chr)"
    return f"{abs(r1p - r2p):,} bp"


def fmt_row(line):
    if not line:
        return ""
    p = line.split("\t")
    if len(p) < 11:
        return line + "\n"
    qname, flag, rname, pos, mapq, cigar = p[0], p[1], p[2], p[3], p[4], p[5]
    rnext, pnext = p[6], p[7] if len(p) > 7 else ""
    seq = p[9][:60] + "..." if len(p[9]) > 60 else p[9]
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


def bhmem_row_is_r1(line: str) -> bool | None:
    p = line.split("\t")
    if len(p) < 2:
        return None
    f = int(p[1])
    if f & 64:
        return True
    if f & 128:
        return False
    return None


def mark_classified_end_bhmem(rows, want_r1: bool):
    out = []
    for line in rows:
        p = line.split("\t")
        if len(p) < 2:
            out.append((line, ""))
            continue
        r1 = bhmem_row_is_r1(line)
        tag = ""
        if r1 is not None and r1 == want_r1:
            tag = "  <<< classified read end (zero-zero compare)"
        out.append((line, tag))
    return out


def mark_classified_end_yap(rows, want_r1: bool):
    out = []
    for line in rows:
        end = parse_yap_read_end(line.split("\t")[0])
        tag = ""
        if end is not None and end == want_r1:
            tag = "  <<< classified read end (zero-zero compare)"
        out.append((line, tag))
    return out


def load_per_read_metrics(path: str) -> dict[tuple[str, str], dict]:
    """(base_id, is_r1 as '0'|'1') -> full row from discrepant per_read TSV."""
    out = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            bid = row["base_id"]
            ir = str(int(row["is_r1"]))
            out[(bid, ir)] = row
    return out


def bhmem_has_both_mates(lines: list[str]) -> bool:
    """Primary R1 and R2 (exclude secondary 0x100)."""
    has_r1 = has_r2 = False
    for line in lines:
        p = line.split("\t")
        if len(p) < 2:
            continue
        f = int(p[1])
        if f & 256:
            continue
        if f & 64:
            has_r1 = True
        if f & 128:
            has_r2 = True
    return has_r1 and has_r2


def yap_has_both_ends(lines: list[str]) -> bool:
    has_r1 = has_r2 = False
    for line in lines:
        end = parse_yap_read_end(line.split("\t")[0])
        if end is True:
            has_r1 = True
        elif end is False:
            has_r2 = True
    return has_r1 and has_r2


def collect_candidates_in_order(
    location_tsv: str,
    metrics: dict[tuple[str, str], dict] | None,
    min_mapq_yap: int,
) -> dict[str, list[dict]]:
    """Per status, ordered list of location rows passing MAPQ (and metrics present)."""
    order = ["same", "same_chrom_diff_span", "diff_chrom", "same_span_diff_strand"]
    pending = {s: [] for s in order}
    with open(location_tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            st = row.get("status", "").strip()
            if st not in pending:
                continue
            bid = row["base_id"]
            ir = str(int(row["is_r1"]))
            if metrics is not None:
                m = metrics.get((bid, ir))
                if not m:
                    continue
                try:
                    mqy = int(m["mapq_yap"])
                except (KeyError, ValueError):
                    continue
                if mqy <= min_mapq_yap:
                    continue
            pending[st].append(row)
    return pending


def pick_examples(
    pending: dict[str, list[dict]],
    per_status: int,
    bhmem_by_base: dict[str, list[str]],
    yap_by_base: dict[str, list[str]],
    require_pair: bool,
) -> dict[str, list[dict]]:
    order = ["same", "same_chrom_diff_span", "diff_chrom", "same_span_diff_strand"]
    buckets = {s: [] for s in order}
    for st in order:
        for row in pending[st]:
            if len(buckets[st]) >= per_status:
                break
            bid = row["base_id"]
            bh = bhmem_by_base.get(bid, [])
            ya = yap_by_base.get(bid, [])
            if require_pair:
                if not bhmem_has_both_mates(bh) or not yap_has_both_ends(ya):
                    continue
            buckets[st].append(row)
    return buckets


def main():
    ap = argparse.ArgumentParser(
        description="PDF: example zero-zero reads by location-match category"
    )
    ap.add_argument(
        "location_compare_tsv",
        help="zero_zero_location_compare.tsv (from compare_zero_zero_locations.py)",
    )
    ap.add_argument("bhmem_bam")
    ap.add_argument("yap_bam", help="yap 3C BAM (Bowtie2)")
    ap.add_argument(
        "-o", "--output", default="zero_zero_example_reads.pdf"
    )
    ap.add_argument(
        "-k",
        "--per-condition",
        type=int,
        default=2,
        help="Number of example read ends per status (default 2)",
    )
    ap.add_argument(
        "--per-read-tsv",
        default="",
        help="discrepant_mismatch_report*.per_read.tsv (MAPQ + mismatch columns). "
        "Recommended: required for yap MAPQ filter and mismatch header.",
    )
    ap.add_argument(
        "--min-mapq-yap",
        type=int,
        default=30,
        help="Keep only read ends with yap MAPQ > this value (default 30).",
    )
    ap.add_argument(
        "--no-require-pair",
        action="store_false",
        dest="require_pair",
        default=True,
        help="Allow single-mate fragments (default: require both R1 and R2 in bhmem and yap).",
    )
    args = ap.parse_args()

    metrics = None
    if args.per_read_tsv:
        metrics = load_per_read_metrics(args.per_read_tsv)
    elif args.min_mapq_yap > 0:
        print(
            "WARNING: no --per-read-tsv; skipping yap MAPQ filter and mismatch header",
            file=sys.stderr,
        )

    pending = collect_candidates_in_order(
        args.location_compare_tsv,
        metrics,
        args.min_mapq_yap,
    )

    base_ids = sorted(
        {
            r["base_id"]
            for lst in pending.values()
            for r in lst
        }
    )
    if not base_ids:
        print("ERROR: no location rows after filters", file=sys.stderr)
        sys.exit(1)

    base_set = set(base_ids)
    yap_by_base: dict[str, list[str]] = defaultdict(list)
    proc = subprocess.Popen(
        ["samtools", "view", args.yap_bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    for line in proc.stdout:
        qname = line.split("\t")[0]
        if "_" not in qname:
            continue
        base = qname.split("_", 1)[0]
        if base in base_set:
            yap_by_base[base].append(line.rstrip())
    proc.wait()

    bhmem_by_base: dict[str, list[str]] = defaultdict(list)
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

    buckets = pick_examples(
        pending,
        args.per_condition,
        bhmem_by_base,
        yap_by_base,
        args.require_pair,
    )

    pages: list[tuple[str, dict]] = []
    for st in ["same", "same_chrom_diff_span", "diff_chrom", "same_span_diff_strand"]:
        for row in buckets[st]:
            pages.append((st, row))

    if not pages:
        print(
            "ERROR: no examples after pair/MAPQ filters; try --no-require-pair or lower --min-mapq-yap",
            file=sys.stderr,
        )
        sys.exit(1)

    status_title = {
        "same": "Same locus (strict): same chr, start, end, strand",
        "same_chrom_diff_span": "Same chromosome, different span",
        "diff_chrom": "Different chromosome (bhmem vs yap)",
        "same_span_diff_strand": "Same span, opposite strand flag",
    }

    with PdfPages(args.output) as pdf:
        for page_i, (status, rec) in enumerate(pages):
            bid = rec["base_id"]
            is_r1 = bool(int(rec["is_r1"]))
            end_label = "R1" if is_r1 else "R2"

            bhmem_rows = bhmem_by_base.get(bid, [])
            yap_rows = yap_by_base.get(bid, [])

            b_r = get_r1_r2_info(bhmem_rows, "bhmem")
            y_r = get_r1_r2_info(yap_rows, "yap")

            blocks = []
            blocks.append(
                f"Page {page_i + 1}/{len(pages)}  |  {status_title.get(status, status)}\n"
            )
            blocks.append(f"Fragment base_id: {bid}  |  Classified read end: {end_label}\n")
            ir_key = str(int(is_r1))
            if metrics:
                m = metrics.get((bid, ir_key))
                if m:
                    blocks.append(
                        "Discrepant per_read (this read end): "
                        f"MAPQ_bhmem={m.get('mapq_bhmem', '')}  MAPQ_yap={m.get('mapq_yap', '')}  "
                        f"NM_bhmem={m.get('nm_bhmem', '')}  "
                        f"yap_raw_mismatch={m.get('yap_raw_mismatch', '')}  "
                        f"yap_corrected_mismatch={m.get('yap_corrected_mismatch', '')}  "
                        f"yap_bisulfite_ignored={m.get('yap_bisulfite_ignored', '')}\n"
                    )
                    blocks.append(
                        "Note: MAPQ and mismatch counts are from the discrepant report (best primary "
                        "alignment per read end). Extra yap rows may be split fragments (-l/-r/-m) "
                        "with lower MAPQ.\n"
                    )
            blocks.append(
                f"TSV: bhmem {rec.get('bhmem_ref', '')}:{rec.get('bhmem_start', '')}-"
                f"{rec.get('bhmem_end', '')} rev={rec.get('bhmem_rev', '')}  |  "
                f"yap {rec.get('yap_ref', '')}:{rec.get('yap_start', '')}-"
                f"{rec.get('yap_end', '')} rev={rec.get('yap_rev', '')}\n"
            )
            if rec.get("max_abs_delta_bp", "") != "":
                blocks.append(
                    f"Δstart {rec.get('delta_start_yap_minus_bhmem', '')}  "
                    f"Δend {rec.get('delta_end_yap_minus_bhmem', '')}  "
                    f"max|Δ| {rec.get('max_abs_delta_bp', '')} bp\n"
                )
            blocks.append("--- Summary: R1-R2 distance (best primary per end, yap) ---\n")
            blocks.append(f"  bhmem: {dist_str(*b_r)}\n")
            blocks.append(f"  yap:   {dist_str(*y_r)}\n")
            blocks.append("=== bhmem BAM ===\n")
            blocks.append(f"Rows: {len(bhmem_rows)}\n")
            for j, (row, tag) in enumerate(
                mark_classified_end_bhmem(bhmem_rows, is_r1)
            ):
                blocks.append(f"--- Row {j+1} ---{tag}\n")
                blocks.append(fmt_row(row))
            blocks.append("\n=== yap 3C BAM ===\n")
            blocks.append(f"Rows: {len(yap_rows)} (full + split with -l/-r/-m)\n")
            for j, (row, tag) in enumerate(
                mark_classified_end_yap(yap_rows, is_r1)
            ):
                blocks.append(f"--- Row {j+1} ---{tag}\n")
                blocks.append(fmt_row(row))
            blocks.append("\n(Split reads: -l left, -r right, -m middle from m3c-split-reads)\n")

            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis("off")
            ax.text(
                0.02,
                0.98,
                "".join(blocks),
                transform=ax.transAxes,
                fontsize=7,
                verticalalignment="top",
                fontfamily="monospace",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.28),
            )
            fig.subplots_adjust(left=0.04, right=0.98, top=0.98, bottom=0.02)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    empty = [s for s in buckets if not buckets[s]]
    print(
        f"Wrote {args.output} ({len(pages)} pages, {args.per_condition} per condition where available)"
    )
    if empty:
        print(f"Note: no examples for: {', '.join(empty)}", file=sys.stderr)


if __name__ == "__main__":
    main()
