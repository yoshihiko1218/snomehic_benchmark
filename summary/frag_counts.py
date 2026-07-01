#!/usr/bin/env python
"""
Uniform fragment-level alignment QC, identical for every dataset.

A FRAGMENT = a distinct read-name (template = one sequenced molecule). For methylation
methods R1 and R2 are aligned independently (possibly in separate BAMs); collapsing by
read-name counts each molecule once (no R1+R2 double-count).

Counts, EXCLUDING duplicates, over one or more per-cell BAMs:
  uniqmap_frag = # distinct read-names with a primary, mapped, non-duplicate alignment
  mapq30_frag  = # of those whose qualifying alignment has MAPQ >= 30
  rate         = 100 * mapq30_frag / uniqmap_frag   (MapQ30 / mapped, fragment-level)

A read qualifies if: not unmapped, not secondary, not supplementary, not duplicate.
(rmdup BAMs already have dups removed; markdup BAMs have them flagged -> is_duplicate.)
Memory-light: streams, keeping a per-read-name best state in a dict only while needed is
avoided by using two sets of read-names. For very large BAMs run per cell (one at a time).

Usage: frag_counts.py --cell <id> --bam a.bam [b.bam ...] [--mapq 30]
Prints: cell  uniqmap_frag  mapq30_frag  rate
"""
import argparse, pysam


def count(bams, mapq=30):
    # BEFORE dedup (include dups) and AFTER dedup (exclude dups), fragment-level.
    uniq_all, hi_all = set(), set()
    uniq_nd, hi_nd = set(), set()
    for bam in bams:
        with pysam.AlignmentFile(bam, "rb") as f:
            for r in f.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                qn = r.query_name
                hq = r.mapping_quality >= mapq
                uniq_all.add(qn)
                if hq:
                    hi_all.add(qn)
                if not r.is_duplicate:
                    uniq_nd.add(qn)
                    if hq:
                        hi_nd.add(qn)
    return len(uniq_all), len(hi_all), len(uniq_nd), len(hi_nd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--bam", nargs="+", required=True)
    ap.add_argument("--mapq", type=int, default=30)
    a = ap.parse_args()
    ua, ha, un, hn = count(a.bam, a.mapq)
    r_all = (100.0 * ha / ua) if ua else float("nan")
    r_nd = (100.0 * hn / un) if un else float("nan")
    # cell  uniq_preDedup  mapq30_preDedup  rate_preDedup  uniq_postDedup  mapq30_postDedup  rate_postDedup
    print(f"{a.cell}\t{ua}\t{ha}\t{r_all:.3f}\t{un}\t{hn}\t{r_nd:.3f}")


if __name__ == "__main__":
    main()
