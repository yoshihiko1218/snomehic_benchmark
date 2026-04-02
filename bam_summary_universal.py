#!/usr/bin/env python3
import argparse
import os
import sys
import pysam


def pct(num, den):
    return 0.0 if den == 0 else (num * 100.0 / den)


def detect_layout(bam_path, check_n=200000):
    """
    Detect whether BAM looks paired-end or single-end by sampling records.
    Returns: "PE", "SE", or "UNKNOWN"
    """
    paired = 0
    unpaired = 0
    checked = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            checked += 1
            if r.is_paired:
                paired += 1
            else:
                unpaired += 1
            if checked >= check_n:
                break

    if checked == 0:
        return "UNKNOWN"

    if paired > 0 and unpaired == 0:
        return "PE"
    if unpaired > 0 and paired == 0:
        return "SE"

    # mixed BAM is unusual; return UNKNOWN so we do not silently assume
    return "UNKNOWN"


def detect_nh_tag(bam_path, check_n=200000):
    """
    Check whether NH tag is present in at least some records.
    """
    checked = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            checked += 1
            if r.has_tag("NH"):
                return True
            if checked >= check_n:
                break
    return False


def should_count_record(r, layout):
    """
    Decide whether this BAM record should represent one counted unit.

    SE: count every record
    PE: count only read1, so counts are fragment-based
    """
    if layout == "SE":
        return True
    elif layout == "PE":
        return r.is_read1
    else:
        return False


def is_primary_alignment(r):
    return (not r.is_unmapped) and (not r.is_secondary) and (not r.is_supplementary)


def summarize_bam(in_bam, out_summary, mapq=30, lambda_name="lambda"):
    if not os.path.exists(in_bam):
        raise FileNotFoundError(f"Input not found: {in_bam}")

    layout = detect_layout(in_bam)
    has_nh = detect_nh_tag(in_bam)

    notes = []

    if layout == "UNKNOWN":
        notes.append(
            "WARNING: could not confidently determine BAM layout (SE vs PE), "
            "or BAM appears mixed. No counts were produced to avoid wrong assumptions."
        )
        with open(out_summary, "w") as out:
            out.write(f"Input:\t{in_bam}\n")
            out.write("Layout:\tUNKNOWN\n")
            out.write(f"NHtagDetected:\t{has_nh}\n\n")
            for n in notes:
                out.write(n + "\n")
        return

    total_units = 0
    primary_mapped = 0
    uniq_mapped = 0
    uniq_mapped_mapq = 0
    uniq_mapped_mapq_to_lambda = 0
    uniq_mapped_mapq_to_target = 0
    duplicate_primary = 0
    secondary_skipped = 0
    supplementary_skipped = 0
    unmapped_skipped = 0
    not_counted_due_to_layout = 0

    with pysam.AlignmentFile(in_bam, "rb") as bam:
        try:
            n_sq = len(bam.header.get("SQ", []))
        except Exception:
            n_sq = 0

        for r in bam.fetch(until_eof=True):
            if not should_count_record(r, layout):
                not_counted_due_to_layout += 1
                continue

            total_units += 1

            if r.is_secondary:
                secondary_skipped += 1
                continue
            if r.is_supplementary:
                supplementary_skipped += 1
                continue
            if r.is_unmapped:
                unmapped_skipped += 1
                continue

            primary_mapped += 1

            if r.is_duplicate:
                duplicate_primary += 1

            # Strict uniqueness only if NH tag exists
            if has_nh:
                try:
                    is_unique = (r.get_tag("NH") == 1)
                except KeyError:
                    is_unique = False
            else:
                # Fallback approximation
                is_unique = True

            if is_unique:
                uniq_mapped += 1

                if r.mapping_quality >= mapq:
                    uniq_mapped_mapq += 1
                    ref = r.reference_name
                    if ref == lambda_name:
                        uniq_mapped_mapq_to_lambda += 1
                    else:
                        uniq_mapped_mapq_to_target += 1

    if total_units == 0:
        notes.append("WARNING: no countable records found.")
    if n_sq == 0:
        notes.append("WARNING: BAM header has no SQ entries.")
    if primary_mapped == 0 and total_units > 0:
        notes.append("WARNING: 0 primary mapped units found.")
    if uniq_mapped_mapq == 0 and uniq_mapped > 0:
        notes.append(f"NOTE: 0 unique units with MAPQ>={mapq}.")
    if not has_nh:
        notes.append(
            "WARNING: NH tag not detected. 'UniqMapped' is approximated as primary mapped units, "
            "not strictly proven unique."
        )
    if layout == "PE":
        notes.append(
            "INFO: Paired-end BAM detected. Counts are fragment-based using read1 only."
        )
    elif layout == "SE":
        notes.append(
            "INFO: Single-end BAM detected. Counts are read-based."
        )

    label_total = "TotalFragments" if layout == "PE" else "TotalReads"
    label_primary = "PrimaryMappedFragments" if layout == "PE" else "PrimaryMappedReads"
    label_uniq = "UniqMappedFragments" if layout == "PE" else "UniqMappedReads"
    label_uniq_mapq = f"UniqMappedMapQ{mapq}Fragments" if layout == "PE" else f"UniqMappedMapQ{mapq}Reads"

    with open(out_summary, "w") as out:
        out.write(f"Input:\t{in_bam}\n")
        out.write(f"Layout:\t{layout}\n")
        out.write(f"NHtagDetected:\t{has_nh}\n")
        out.write(f"{label_total}:\t{total_units}\t100%\n")
        out.write(f"{label_primary}:\t{primary_mapped}\t{pct(primary_mapped, total_units):.6f}%\n")
        out.write(f"{label_uniq}:\t{uniq_mapped}\t{pct(uniq_mapped, total_units):.6f}%\n")
        out.write(f"{label_uniq_mapq}:\t{uniq_mapped_mapq}\t{pct(uniq_mapped_mapq, total_units):.6f}%\n")
        out.write(
            f"{label_uniq_mapq}ToLambda:\t{uniq_mapped_mapq_to_lambda}\t"
            f"{pct(uniq_mapped_mapq_to_lambda, uniq_mapped_mapq):.6f}%\n"
        )
        out.write(
            f"{label_uniq_mapq}ToTargetSpecies:\t{uniq_mapped_mapq_to_target}\t"
            f"{pct(uniq_mapped_mapq_to_target, uniq_mapped_mapq):.6f}%\n"
        )
        out.write(
            f"DuplicatePrimaryMapped:\t{duplicate_primary}\t"
            f"{pct(duplicate_primary, primary_mapped):.6f}%\n"
        )
        out.write(
            f"SecondarySkipped:\t{secondary_skipped}\t"
            f"{pct(secondary_skipped, total_units):.6f}%\n"
        )
        out.write(
            f"SupplementarySkipped:\t{supplementary_skipped}\t"
            f"{pct(supplementary_skipped, total_units):.6f}%\n"
        )
        out.write(
            f"UnmappedSkipped:\t{unmapped_skipped}\t"
            f"{pct(unmapped_skipped, total_units):.6f}%\n"
        )
        out.write(f"NotCountedDueToLayout:\t{not_counted_due_to_layout}\n")

        if notes:
            out.write("\n")
            for n in notes:
                out.write(n + "\n")


def main():
    p = argparse.ArgumentParser(
        description="Universal BAM QC summary for SE/PE BAMs with conservative sanity checks."
    )
    p.add_argument("--in_bam", required=True, help="Input BAM")
    p.add_argument("--out_summary", required=True, help="Output summary txt")
    p.add_argument("--mapq", type=int, default=30, help="MAPQ cutoff (default: 30)")
    p.add_argument(
        "--lambda_name",
        default="lambda",
        help='Reference name for lambda spike-in (default: "lambda")'
    )
    args = p.parse_args()

    try:
        summarize_bam(
            in_bam=args.in_bam,
            out_summary=args.out_summary,
            mapq=args.mapq,
            lambda_name=args.lambda_name,
        )
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
