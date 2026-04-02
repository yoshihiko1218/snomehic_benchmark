
#!/usr/bin/env python3
import argparse
import os
import pysam

def pct(num, den):
    return 0.0 if den == 0 else (num * 100.0 / den)

def se_bam_summary(in_bam, out_summary, mapq=30, lambda_name="lambda"):
    if not os.path.exists(in_bam):
        raise FileNotFoundError(f"Input not found: {in_bam}")

    total = 0
    mapped = 0
    mapped_mapq = 0
    mapped_mapq_to_lambda = 0
    mapped_mapq_to_target = 0
    dup_mapped = 0
    secondary = 0
    supplementary = 0

    with pysam.AlignmentFile(in_bam, "rb") as bam:
        # sanity: has reference sequences in header?
        try:
            n_sq = len(bam.header.get("SQ", []))
        except Exception:
            n_sq = 0

        for r in bam.fetch(until_eof=True):
            total += 1

            if r.is_secondary:
                secondary += 1
                continue
            if r.is_supplementary:
                supplementary += 1
                continue

            if r.is_unmapped:
                continue

            mapped += 1

            if r.is_duplicate:
                dup_mapped += 1

            if r.mapping_quality >= mapq:
                mapped_mapq += 1
                ref = r.reference_name
                if ref == lambda_name:
                    mapped_mapq_to_lambda += 1
                else:
                    mapped_mapq_to_target += 1

    # sanity: warn-style notes in output if something looks off
    notes = []
    if total == 0:
        notes.append("WARNING: no alignments iterated (empty file or unreadable).")
    if n_sq == 0:
        notes.append("WARNING: BAM header has no SQ entries (unexpected).")
    if mapped == 0 and total > 0:
        notes.append("WARNING: 0 mapped reads found; check mapping or flags.")
    if mapped_mapq == 0 and mapped > 0:
        notes.append(f"NOTE: 0 reads with MAPQ>={mapq}; might be expected for this aligner/settings.")

    with open(out_summary, "w") as out:
        out.write(f"Input:\t{in_bam}\n")
        out.write(f"TotalReads:\t{total}\t100%\n")
        out.write(f"Mapped:\t{mapped}\t{pct(mapped, total):.4f}%\n")
        out.write(f"MappedMapQ{mapq}:\t{mapped_mapq}\t{pct(mapped_mapq, total):.4f}%\n")
        out.write(f"MappedMapQ{mapq}ToLambda:\t{mapped_mapq_to_lambda}\t{pct(mapped_mapq_to_lambda, mapped_mapq):.4f}%\n")
        out.write(f"MappedMapQ{mapq}ToTargetSpecies:\t{mapped_mapq_to_target}\t{pct(mapped_mapq_to_target, mapped_mapq):.4f}%\n")
        out.write(f"DuplicateMapped:\t{dup_mapped}\t{pct(dup_mapped, mapped):.4f}%\n")
        out.write(f"SecondarySkipped:\t{secondary}\t{pct(secondary, total):.4f}%\n")
        out.write(f"SupplementarySkipped:\t{supplementary}\t{pct(supplementary, total):.4f}%\n")
        if notes:
            out.write("\n")
            for n in notes:
                out.write(n + "\n")

def main():
    p = argparse.ArgumentParser(description="Single-end BAM summary (non-Hi-C metrics).")
    p.add_argument("--in_bam", required=True, help="Input BAM (single-end ok).")
    p.add_argument("--out_summary", required=True, help="Output summary txt.")
    p.add_argument("--mapq", type=int, default=30, help="MAPQ cutoff (default: 30).")
    p.add_argument("--lambda_name", default="lambda", help='Reference name for lambda (default: "lambda").')
    args = p.parse_args()
    se_bam_summary(args.in_bam, args.out_summary, mapq=args.mapq, lambda_name=args.lambda_name)

if __name__ == "__main__":
    main()

