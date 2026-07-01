#!/usr/bin/env python
"""
snmCAT has no MapQ30 columns in its YAP MappingSummary, so compute the MapQ30 rate
directly from the correct BAM = the DNA-reads BAM (mapping_brain/Group*/bam/<cell>.dna_reads.bam).
For each cell, count primary alignments (-F 0x904) total and with MAPQ>=30, and pull
R1/R2 InputReads + UniqueMappedReads from MappingSummary so the denominator can match
the other datasets' "MapQ30 / input" definition.
Writes summary/snmCAT_mapq30_percell.tsv:
  cell  mapq30_n  primary_mapped_n  input_reads  uniq_reads  mapq30_rate_vs_input  mapq30_rate_vs_mapped
"""
import glob, os, subprocess
import pandas as pd

BASE = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
bams = sorted(glob.glob(f"{BASE}/snmCAT/mapping_brain/Group*/bam/*.dna_reads.bam"))
ms = pd.read_csv(f"{BASE}/snmCAT/mapping_brain/stats/MappingSummary.csv.gz",
                 compression="gzip", index_col=0)


def count(bam, q=None):
    cmd = ["samtools", "view", "-c", "-F", "0x904", bam]
    if q is not None:
        cmd[3:3] = ["-q", str(q)]
    return int(subprocess.run(cmd, capture_output=True, text=True).stdout.strip() or 0)


rows = []
for i, bam in enumerate(bams, 1):
    cell = os.path.basename(bam)[: -len(".dna_reads.bam")]
    mapped = count(bam)
    mq30 = count(bam, 30)
    inp = uniq = float("nan")
    if cell in ms.index:
        r = ms.loc[cell]
        inp = float(r.get("R1InputReads", 0)) + float(r.get("R2InputReads", 0))
        uniq = float(r.get("R1UniqueMappedReads", 0)) + float(r.get("R2UniqueMappedReads", 0))
    rows.append({
        "cell": cell, "mapq30_n": mq30, "primary_mapped_n": mapped,
        "input_reads": inp, "uniq_reads": uniq,
        "mapq30_rate_vs_input": (100.0 * mq30 / inp) if inp and inp == inp and inp > 0 else float("nan"),
        "mapq30_rate_vs_mapped": (100.0 * mq30 / mapped) if mapped > 0 else float("nan"),
    })
    if i % 20 == 0:
        print(f"  {i}/{len(bams)}")

df = pd.DataFrame(rows)
out = f"{BASE}/summary/snmCAT_mapq30_percell.tsv"
df.to_csv(out, sep="\t", index=False)
print("wrote", out, len(df), "cells")
print(df[["mapq30_n", "primary_mapped_n", "input_reads", "uniq_reads",
          "mapq30_rate_vs_input", "mapq30_rate_vs_mapped"]].describe().round(2).to_string())
