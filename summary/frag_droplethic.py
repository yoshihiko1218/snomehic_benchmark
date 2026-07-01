#!/usr/bin/env python
"""
Per-barcode fragment-level QC for droplet Hi-C, consistent with frag_counts.py:
stream the name-sorted CB-tagged BAM, pair consecutive mates by read-name, and for
each barcode count (EXCLUDING duplicates, OR logic = fragment qualifies if >=1 mate
is primary/mapped/non-dup):
  uniqmap_frag : # fragments with >=1 primary, mapped, non-dup mate
  mapq30_frag  : # of those with >=1 such mate at MAPQ>=30
Then restrict to the valid-barcode set and write summary/frag_counts/droplethic_percell.tsv.
"""
import argparse, collections, pysam, pandas as pd

def get_cb(r):
    try:
        return r.get_tag("CB")
    except KeyError:
        qn = r.query_name
        return qn.rsplit(":", 1)[-1] if ":" in qn else "NA"

def states(r, mapq):
    # returns (mapped_inclDup, mapq30_inclDup, mapped_noDup, mapq30_noDup) for one mate
    mapped = (not r.is_unmapped and not r.is_secondary and not r.is_supplementary)
    hq = mapped and r.mapping_quality >= mapq
    nd = mapped and not r.is_duplicate
    hq_nd = nd and r.mapping_quality >= mapq
    return mapped, hq, nd, hq_nd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--valid", required=True, help="valid per_cell tsv (CB column)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mapq", type=int, default=30)
    a = ap.parse_args()

    uniq_all, hi_all = collections.Counter(), collections.Counter()
    uniq_nd, hi_nd = collections.Counter(), collections.Counter()
    prev = None
    n = 0
    with pysam.AlignmentFile(a.bam, "rb") as f:
        for r in f.fetch(until_eof=True):
            n += 1
            if prev is None:
                prev = r; continue
            if prev.query_name != r.query_name:
                prev = r; continue   # unpaired/orphan; advance
            r1, r2 = prev, r
            prev = None
            cb = get_cb(r1)
            m1, h1, nd1, hnd1 = states(r1, a.mapq)
            m2, h2, nd2, hnd2 = states(r2, a.mapq)
            if m1 or m2:
                uniq_all[cb] += 1
                if h1 or h2:
                    hi_all[cb] += 1
            if nd1 or nd2:
                uniq_nd[cb] += 1
                if hnd1 or hnd2:
                    hi_nd[cb] += 1
            if n % 50_000_000 == 0:
                print(f"  {n:,} alignments", flush=True)

    valid = set(pd.read_csv(a.valid, sep="\t")["CB"].astype(str))
    rows = []
    for cb in valid:
        ua, ha = uniq_all.get(cb, 0), hi_all.get(cb, 0)
        un, hn = uniq_nd.get(cb, 0), hi_nd.get(cb, 0)
        rows.append({"cell": cb,
                     "uniq_preDedup": ua, "mapq30_preDedup": ha,
                     "rate_preDedup": (100.0 * ha / ua) if ua else float("nan"),
                     "uniq_postDedup": un, "mapq30_postDedup": hn,
                     "rate_postDedup": (100.0 * hn / un) if un else float("nan")})
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"wrote {a.out}: {len(rows)} valid barcodes (of {len(uniq_all):,} total)")

if __name__ == "__main__":
    main()
