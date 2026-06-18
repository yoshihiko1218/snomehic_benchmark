import gzip, subprocess

HG = "/gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38_bismark/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
SAM = "/gpfs/projects/b1198/epifluidlab/yoshii/software/samtools-1.16/bin/samtools"
METHY = "/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/scnome/05.methy"
CELL = "K562_01"
CHROM = "chr21"

out = subprocess.run([SAM, "faidx", HG, CHROM], capture_output=True, text=True).stdout
seq = "".join(out.split("\n")[1:]).upper()
print(f"{CHROM} length: {len(seq):,}")

def scan(path):
    n_cpg_cyt = n_gpc_cyt = 0
    hcg, gch = set(), set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t", 2)
            if f[0] != CHROM:
                continue
            i = int(f[1]) - 1
            if i < 2 or i + 2 >= len(seq):
                continue
            b = seq[i]
            # CpG context cytosine?
            if b == "C" and seq[i + 1] == "G":
                cpg = i; n_cpg_cyt += 1
            elif b == "G" and seq[i - 1] == "C":
                cpg = i - 1; n_cpg_cyt += 1
            else:
                cpg = -1
            if cpg >= 1 and seq[cpg - 1] != "G":
                hcg.add((f[0], cpg))
            # GpC context cytosine?
            if b == "C" and seq[i - 1] == "G":
                g = i - 1; is_gcg = seq[i + 1] == "G"; n_gpc_cyt += 1
            elif b == "G" and seq[i + 1] == "C":
                g = i; is_gcg = seq[i + 2] == "G"; n_gpc_cyt += 1
            else:
                g = -1; is_gcg = True
            if g >= 0 and not is_gcg:
                gch.add((f[0], g))
    return n_cpg_cyt, n_gpc_cyt, hcg, gch

c1, g1, hcg1, gch1 = scan(f"{METHY}/{CELL}_1.rmdup.bismark.cov.gz")
c2, g2, hcg2, gch2 = scan(f"{METHY}/{CELL}_2.rmdup.bismark.cov.gz")

print(f"\n=== (B) STRAND DESTRANDING -- +/- of same CpG/GpC merged to ONE locus ({CELL}_1, {CHROM}) ===")
print(f"HCG: {c1:,} covered CpG-context cytosine lines (both strands)  ->  {len(hcg1):,} unique loci"
      f"   => {c1/max(len(hcg1),1):.2f} strand-observations per locus")
print(f"GCH: {g1:,} covered GpC-context cytosine lines (both strands)  ->  {len(gch1):,} unique loci"
      f"   => {g1/max(len(gch1),1):.2f} strand-observations per locus")
print("  (ratio >1 confirms + and - strand cytosines of the same site collapse into one locus)")

print(f"\n=== (A) MATE UNION -- a locus seen in BOTH R1 and R2 counted ONCE ({CHROM}) ===")
for name, m1, m2 in [("HCG", hcg1, hcg2), ("GCH", gch1, gch2)]:
    inter = len(m1 & m2); uni = len(m1 | m2)
    print(f"{name}: |R1|={len(m1):,}  |R2|={len(m2):,}  |overlap|={inter:,}  |UNION|={uni:,}")
    print(f"     naive sum would be {len(m1)+len(m2):,}  (would double-count {inter:,} shared loci)")
    print(f"     union == |R1|+|R2|-|overlap| ? {uni == len(m1)+len(m2)-inter}")
