#!/usr/bin/env python
"""
Build per-dataset manifests for the uniform fragment-level QC recompute.
Each manifest line: cell <TAB> bam1[,bam2,...] <TAB> needs_markdup(0/1)
Only cells whose BAM(s) exist are written; coverage is printed for sanity.
"""
import os, glob, pandas as pd

B = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark"
GM = "/home/jmj7858/epifluidlab/sc_nomehic_cellline/gm_sc_new/04.alignment_snakemake"
OUT = f"{B}/summary/frag_jobs"
os.makedirs(OUT, exist_ok=True)


def write(ds, rows):
    rows = [r for r in rows if all(os.path.exists(p) for p in r[1].split(","))]
    with open(f"{OUT}/{ds}.manifest.tsv", "w") as fh:
        for cell, bams, md in rows:
            fh.write(f"{cell}\t{bams}\t{md}\n")
    print(f"{ds:11s}: {len(rows)} cells -> {ds}.manifest.tsv")


# nagano (15) — markdup.bam, dup-flagged
cells = [c.strip() for c in open(f"{B}/nagano/acc_list.txt") if c.strip()]
write("nagano", [(c, f"{B}/nagano/alignment/{c}.markdup.bam", 0) for c in cells])

# smallwood (51) — rmdup.RG.bam (dup-flagged), SE
cells = pd.read_csv(f"{B}/smallwood/smallwood_qc_summary.csv")["CellID"].astype(str)
write("smallwood", [(c, f"{B}/smallwood/05.align_mm10/{c}.rmdup.RG.bam", 0) for c in cells])

# scnome (23) — per-mate rmdup.RG.bam (dup-flagged)
cells = pd.read_csv(f"{B}/scnome/scnome_qc_summary.csv")["CellID"].astype(str)
write("scnome", [(c, f"{B}/scnome/04.alignment/{c}_1.rmdup.RG.bam,{B}/scnome/04.alignment/{c}_2.rmdup.RG.bam", 0) for c in cells])

# snmCseq3 — calmd.bam (dup-flagged)
bams = sorted(glob.glob(f"{B}/snmCseq3/04.bhmem_bam/*.calmd.bam"))
write("snmCseq3", [(os.path.basename(b)[:-len('.calmd.bam')], b, 0) for b in bams])

# scnomehic (187 passed) — external calmd.bam (dup-flagged)
passed = set(pd.read_csv("/home/jmj7858/epifluidlab/scnomehic_paper/s1/QC/gm_passed.txt", header=None)[0].astype(str))
write("scnomehic", [(c, f"{GM}/{c}.calmd.bam", 0) for c in sorted(passed)])

# snmCseq2 (96 mm10) — raw clean_bismark_bt2.bam per mate -> NEEDS MARKDUP
mm10 = pd.read_csv(f"{B}/snmCseq2/codes/cells_mm10.txt", header=None)[0].astype(str)
write("snmCseq2", [(c, f"{B}/snmCseq2/05.align/{c}_1.clean_bismark_bt2.bam,{B}/snmCseq2/05.align/{c}_2.clean_bismark_bt2.bam", 1) for c in mm10])

# snmCAT — dna_reads.bam -> NEEDS MARKDUP
bams = sorted(glob.glob(f"{B}/snmCAT/mapping_brain/Group*/bam/*.dna_reads.bam"))
write("snmCAT", [(os.path.basename(b)[:-len('.dna_reads.bam')], b, 1) for b in bams])
