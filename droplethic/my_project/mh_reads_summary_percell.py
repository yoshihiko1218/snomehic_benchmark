import pysam
import argparse
from collections import defaultdict

def get_cb(read):
    # Prefer CB tag (what your BAM has)
    if read.has_tag("CB"):
        return read.get_tag("CB")
    # Fallback: parse from QNAME suffix after ":" (matches your example)
    qn = read.query_name
    return qn.split(":")[-1] if ":" in qn else "NO_CB"

def mh_reads_summary_percell(in_bam, out_tsv, mapq=30, lambda_name="lambda"):
    # counts[CB][metric] = int
    counts = defaultdict(lambda: defaultdict(int))

    mode = "rb" if (in_bam == "-" or in_bam.endswith(".bam")) else "rc"
    with pysam.AlignmentFile(in_bam, mode) as samfile:
        pre_read = None
        for read in samfile.fetch(until_eof=True):
            if pre_read is None:
                pre_read = read
                continue

            # enforce adjacent mates by QNAME
            if pre_read.query_name != read.query_name:
                pre_read = read
                continue

            r1, r2 = pre_read, read
            pre_read = None

            cb = get_cb(r1)
            c = counts[cb]

            c["TotalFragments"] += 1

            # original "UniqMapped" definition (primary, paired, both mapped)
            if (r1.is_paired and (not r1.is_unmapped) and (not r2.is_unmapped)
                and (not r1.is_secondary) and (not r2.is_secondary)):
                c["UniqMapped"] += 1

                if r1.mapping_quality >= mapq and r2.mapping_quality >= mapq:
                    c["UniqMappedMapQ30"] += 1

                    # NOTE: original code checks read.reference_name == "lambda" using only one mate;
                    # here we follow same spirit but use either mate (safer).
                    if (r1.reference_name == lambda_name) or (r2.reference_name == lambda_name):
                        c["UniqMappedMapQ30ToLambda"] += 1
                    else:
                        c["UniqMappedMapQ30ToTargetSpecies"] += 1

                        # NoPcr branch: requires duplicate flag to be set in BAM
                        if (not r1.is_duplicate) and (not r2.is_duplicate):
                            c["UniqMappedMapQ30NoPcr"] += 1
                            c["UniqMappedMapQ30NoPcrToTargetSpecies"] += 1

                            # cis/trans like original (based on reference vs next_reference)
                            # Use r1 vs its mate reference; pysam sets next_reference_name.
                            # This is equivalent to your original check.
                            if r1.reference_name != r1.next_reference_name:
                                c["UniqMappedMapQ30NoPcrTrans"] += 1
                            else:
                                c["UniqMappedMapQ30NoPcrCis"] += 1
                                tlen = abs(r1.template_length)
                                if tlen >= 1000:
                                    c["UniqMappedMapQ30NoPcrCisMore1kb"] += 1
                                    if tlen >= 20000:
                                        c["UniqMappedMapQ30NoPcrCisMore20kb"] += 1
                # else: mapped but low MAPQ -> not counted further

    # Write per-cell table (TSV)
    header = [
        "CB",
        "TotalFragments",
        "UniqMapped",
        "UniqMappedMapQ30",
        "UniqMappedMapQ30ToLambda",
        "UniqMappedMapQ30ToTargetSpecies",
        "UniqMappedMapQ30NoPcr",
        "UniqMappedMapQ30NoPcrToLambda",
        "UniqMappedMapQ30NoPcrToTargetSpecies",
        "UniqMappedMapQ30NoPcrTrans",
        "UniqMappedMapQ30NoPcrCis",
        "UniqMappedMapQ30NoPcrCisMore1kb",
        "UniqMappedMapQ30NoPcrCisMore20kb",
        # Useful ratios
        "MapQ30_over_TotalFragments",
        "MapQ30_over_UniqMapped",
        "TransFrac_over_TargetNoPcr",
        "CisFrac_over_TargetNoPcr",
    ]

    with open(out_tsv, "w") as out:
        out.write("\t".join(header) + "\n")
        for cb in sorted(counts.keys()):
            c = counts[cb]
            total = c.get("TotalFragments", 0)
            um = c.get("UniqMapped", 0)
            mq = c.get("UniqMappedMapQ30", 0)
            targ_nopcr = c.get("UniqMappedMapQ30NoPcrToTargetSpecies", 0)
            trans = c.get("UniqMappedMapQ30NoPcrTrans", 0)
            cis = c.get("UniqMappedMapQ30NoPcrCis", 0)

            mapq_over_total = (mq / total) if total else 0.0
            mapq_over_um = (mq / um) if um else 0.0
            trans_frac = (trans / targ_nopcr) if targ_nopcr else 0.0
            cis_frac = (cis / targ_nopcr) if targ_nopcr else 0.0

            row = [
                cb,
                str(total),
                str(um),
                str(mq),
                str(c.get("UniqMappedMapQ30ToLambda", 0)),
                str(c.get("UniqMappedMapQ30ToTargetSpecies", 0)),
                str(c.get("UniqMappedMapQ30NoPcr", 0)),
                str(c.get("UniqMappedMapQ30NoPcrToLambda", 0)),
                str(targ_nopcr),
                str(trans),
                str(cis),
                str(c.get("UniqMappedMapQ30NoPcrCisMore1kb", 0)),
                str(c.get("UniqMappedMapQ30NoPcrCisMore20kb", 0)),
                f"{mapq_over_total:.6f}",
                f"{mapq_over_um:.6f}",
                f"{trans_frac:.6f}",
                f"{cis_frac:.6f}",
            ]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--in_bam", required=True, help="BAM/CRAM path or '-' for stdin. MUST be read-name sorted stream.")
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--mapq", type=int, default=30)
    p.add_argument("--lambda_name", default="lambda")
    args = p.parse_args()
    mh_reads_summary_percell(args.in_bam, args.out_tsv, mapq=args.mapq, lambda_name=args.lambda_name)
