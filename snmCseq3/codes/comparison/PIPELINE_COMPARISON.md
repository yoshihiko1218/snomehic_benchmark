# Bhmem vs Bismark (Yap) Pipeline Comparison for snmCseq3

Detailed step-by-step comparison of the two alignment pipelines as actually run in this benchmark, with concrete examples showing how read conversion and genome alignment work internally.

## Table of Contents

- [Actual Flags Used](#actual-flags-used)
- [Bismark (Yap Pipeline)](#bismark-yap-pipeline)
  - [R1: bismark --pbat](#r1-bismark---pbat-single-end)
  - [R2: bismark standard directional](#r2-bismark-standard-directional-single-end)
- [Bhmem Pipeline](#bhmem-pipeline)
  - [Read Conversion](#step-2a-read-conversion)
  - [Four Alignment Trials](#step-2b-four-alignment-trials)
  - [SAM Record Construction](#step-2c-sam-record-construction-string2samrecord)
  - [Pairing with -nonDirectional](#step-2d-pairing-with--nondirectional)
  - [Best Pair Selection](#step-2e-best-pair-selection)
- [Side-by-Side Comparison with Concrete Example](#side-by-side-comparison-with-concrete-example)
- [Summary of Equivalences and Differences](#summary-of-equivalences-and-differences)

---

## Actual Flags Used

### Bhmem (`02.alignment.sh`)

```bash
java -Xmx60G ... Bhmem mm10.fa out.bam R1.fq.gz R2.fq.gz \
  -rgId {prefix} -rgSm snm3C \
  -nonDirectional -pbat \
  -buffer 100000 \
  -enzymeList dpnII.span_region.bedgraph \
  -outputMateDiffChr
```

Key flags: **`-pbat`** and **`-nonDirectional`** (all 4 cross-pair combinations tried).

### Bismark (yap Snakefile)

```bash
# R1: PBAT mode
bismark {bismark_reference} -un --bowtie2 {input} --pbat -o bam/ --temp_dir bam/

# R2: Standard directional mode (no --pbat)
bismark {bismark_reference} -un --bowtie2 {input} -o bam/ --temp_dir bam/
```

Key: **R1 uses `--pbat`**, **R2 uses standard directional (no --pbat)**.

---

## Bismark (Yap Pipeline)

Source: `/home/jmj7858/epifluidlab/software/Bismark-0.24.2/bismark`

Bismark maps R1 and R2 as **independent single-end** reads with different strand modes.

### Genome Indices

Bismark's `bismark_genome_preparation` creates two converted genomes:
- **CT genome**: every C in reference replaced with T (`Bisulfite_Genome/CT_conversion/`)
- **GA genome**: every G in reference replaced with A (`Bisulfite_Genome/GA_conversion/`)

Both are indexed by Bowtie2.

### R1: `bismark --pbat` (single-end)

#### Step 1: Read conversion

For single-end `--pbat`, bismark creates **only a G->A converted** query file (bismark line 528-535):

```perl
elsif($pbat){
    ($G_to_A_infile) = biTransformFastQFiles ($filename);
    $fhs[0]->{inputfile} = $fhs[1]->{inputfile} = $G_to_A_infile;
}
```

The conversion is a simple global substitution (line ~5495):

```perl
$sequence_G_to_A = $sequence;
$sequence_G_to_A =~ tr/G/A/;
```

**Example:**
```
Original R1:  ACGTCGATCG
G->A convert: ACATCAATCA
```

#### Step 2: Genome index initialization

For single-end `--pbat`, only 2 filehandles are created (bismark line 7171-7185):

```perl
@fhs = (
  { name => 'GAreadCTgenome', bisulfiteIndex => $CT_index_basename },  # index 0
  { name => 'GAreadGAgenome', bisulfiteIndex => $GA_index_basename },  # index 1
);
```

Both receive the **same** G->A converted query file as input.

#### Step 3: Two alignment trials

Bowtie2 is launched twice with orientation restrictions (bismark line 6845-6850):
- `GAreadCTgenome` gets `--nofw` (only align to reverse complement of reference)
- `GAreadGAgenome` gets `--norc` (only align forward to reference)

| Trial | Converted query | Genome index | Bowtie2 flag | Bismark strand |
|-------|----------------|-------------|-------------|---------------|
| 1 | `ACATCAATCA` (G->A) | CT genome | `--nofw` | CTOT |
| 2 | `ACATCAATCA` (G->A) | GA genome | `--norc` | CTOB |

#### Step 4: Best alignment selection

From the 2 trials, bismark collects all alignments keyed by `chr:position:strand` (bismark line 3004-3060):

- If exactly **1 unique position** found -> accept it
- If **2+ positions** with **different** AS scores -> keep the one with **highest AS**
- If **2+ positions** with **tied** AS scores -> **reject as ambiguous**

#### Step 5: Index remapping for methylation calling

Because `--pbat` is set, a modifier shifts the internal index (bismark line 4280-4284):

```perl
if ($pbat){
    $pbat_index_modifier += 2;
}
```

The @fhs indices 0,1 become conceptual indices 2,3:
- Index 2 (CTOT): `XR:Z:GA`, `XG:Z:CT`, strand = `-`
- Index 3 (CTOB): `XR:Z:GA`, `XG:Z:GA`, strand = `+`

---

### R2: `bismark` standard directional (single-end)

#### Step 1: Read conversion

Standard directional mode creates **only a C->T converted** query file (bismark line 519-526):

```perl
if ($directional){
    ($C_to_T_infile) = biTransformFastQFiles ($filename);
    $fhs[0]->{inputfile} = $fhs[1]->{inputfile} = $C_to_T_infile;
}
```

**Example:**
```
Original R2:  TGCATGCATG
C->T convert: TGTATGTATG
```

#### Step 2: Genome index initialization

For single-end directional, 2 filehandles (bismark line 7125-7139):

```perl
@fhs = (
  { name => 'CTreadCTgenome', bisulfiteIndex => $CT_index_basename },  # index 0
  { name => 'CTreadGAgenome', bisulfiteIndex => $GA_index_basename },  # index 1
);
```

#### Step 3: Two alignment trials

| Trial | Converted query | Genome index | Bowtie2 flag | Bismark strand |
|-------|----------------|-------------|-------------|---------------|
| 1 | `TGTATGTATG` (C->T) | CT genome | `--norc` (forward only) | OT |
| 2 | `TGTATGTATG` (C->T) | GA genome | `--nofw` (revcomp only) | OB |

#### Step 4: Best alignment selection

Same as R1: unique best by AS, reject ties.

#### Step 5: Index determines tags (no pbat modifier)

- Index 0 (OT): `XR:Z:CT`, `XG:Z:CT`, strand = `+`
- Index 1 (OB): `XR:Z:CT`, `XG:Z:GA`, strand = `-`

---

### Post-alignment steps (yap Snakefile)

1. **Split-read remapping**: Unmapped reads are split into left (40bp), middle (>=30bp), right (40bp) fragments by `yap-internal m3c-split-reads`, then remapped with the same bismark flags per mate
2. **Merge**: first-pass BAM + split-read BAM per mate
3. **MAPQ filter**: `samtools view -q 10`
4. **Sort + dedup**: coordinate sort, then `picard MarkDuplicates REMOVE_DUPLICATES=true`
5. **3C BAM** (used for comparison): merge pre-dedup, post-MAPQ-filter R1+R2, name-sort -> `{cell_id}.3C.sorted.bam`

---

## Bhmem Pipeline

Source: `/home/jmj7858/epifluidlab/software/bisulfitehic/src/main/java/edu/mit/compbio/bisulfitehic/mapping/Bhmem.java`

Bhmem maps R1 and R2 as **paired-end** with internal BWA mem alignment.

### Step 1: Trimming (`01.trim.sh`)

Same cutadapt parameters as yap:
- Pass 1: adapter removal
- Pass 2: `-O 6 -q 20 -u 10 -u -10 -m 30` (trim 10bp from both ends, Q>=20, min length 30)

**No read renaming** (unlike yap which appends `_1`/`_2` suffix).

### Step 2a: Read conversion

(Bhmem.java line 281-290)

```java
boolean pbatR1 = snm3c ? true : pbat;   // pbat=true, snm3c=false -> pbatR1=true
boolean pbatR2 = snm3c ? false : pbat;  // pbat=true, snm3c=false -> pbatR2=true

// pbatR1=true -> R1 gets G->A
ShortRead read1ct = new ShortRead(name, seq1.replace('G','A').getBytes(), quals);

// pbatR2=true -> R2 gets C->T
ShortRead read2ga = new ShortRead(name, seq2.replace('C','T').getBytes(), quals);
```

**Example:**
```
Original R1:  ACGTCGATCG
R1 query:     ACATCAATCA   (G->A)

Original R2:  TGCATGCATG
R2 query:     TGTATGTATG   (C->T)
```

Note: the variable names `read1ct` and `read2ga` are misleading -- they are inherited from non-PBAT mode and don't reflect what conversion was actually applied.

The **original unconverted sequences** are stored in lists `L1` and `L2` for later restoration.

### Step 2b: Four alignment trials

(Bhmem.java line 251-258)

```java
// pbatR1=true -> genome indices are SWAPPED for R1
String[] samsEnd1CT = memGA.align(L1CT);  // R1 G->A query -> GA genome
String[] samsEnd1GA = memCT.align(L1CT);  // R1 G->A query -> CT genome

// pbatR2=true -> genome indices are SWAPPED for R2
String[] samsEnd2CT = memGA.align(L2GA);  // R2 C->T query -> GA genome
String[] samsEnd2GA = memCT.align(L2GA);  // R2 C->T query -> CT genome
```

| Trial variable | Mate | Converted query | Genome index |
|----------------|------|----------------|-------------|
| `samsEnd1CT` | R1 | `ACATCAATCA` (G->A) | **GA** genome |
| `samsEnd1GA` | R1 | `ACATCAATCA` (G->A) | **CT** genome |
| `samsEnd2CT` | R2 | `TGTATGTATG` (C->T) | **GA** genome |
| `samsEnd2GA` | R2 | `TGTATGTATG` (C->T) | **CT** genome |

Note: the trial variable names (`End1CT`, `End1GA`) are also misleading due to the PBAT swap. `samsEnd1CT` actually aligns to the **GA** genome when `pbat=true`.

Each trial runs BWA mem in single-end mode. BWA can produce multiple hits (primary + secondary/supplementary). BWA computes MAPQ, AS, NM, and CIGAR against the converted reference.

### Step 2c: SAM record construction (`String2SamRecord`)

(Bhmem.java line 620-692)

For each BWA SAM output line, bhmem builds a SAMRecord through these transformations:

#### Example: R1 mapped forward (flag=0) to GA genome at chr1:1000

BWA output:
```
read1  0  chr1_GA_converted  1000  60  10M  *  0  0  ACATCAATCA  IIIIIIIIII  NM:i:0  AS:i:0
```

**a) Strip genome suffix** (line 631-632):
```java
splitin[2] = splitin[2].replace("_CT_converted", "");
splitin[2] = splitin[2].replace("_GA_converted", "");
// chr1_GA_converted -> chr1
```

**b) PBAT strand flip** (line 642-652, `pbatThisEnd=true` for R1):
```java
if (pbatThisEnd) {
    record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());  // false -> true
    originalSeq = SequenceUtil.reverseComplement(originalSeq);  // ACGTCGATCG -> CGATCGACGT
    baseQ = new StringBuilder(baseQ).reverse().toString();  // reversed
}
```

**c) Standard BAM reverse-complement** (line 655-658, now `readNegativeStrand=true`):
```java
if (record.getReadNegativeStrandFlag()) {
    originalSeq = SequenceUtil.reverseComplement(originalSeq);  // CGATCGACGT -> ACGTCGATCG
    baseQ = new StringBuilder(baseQ).reverse().toString();  // reversed back
}
```

**d) Set sequence from original** (line 659):
```java
record.setReadString(modifySeqByCigar(originalSeq, record.getCigarString()));
// SEQ = ACGTCGATCG (original restored after double revcomp)
```

**Net effect for R1 mapped forward by BWA:**
- Double revcomp cancels: SEQ = original sequence
- Strand flag flipped: `false -> true` (marked reverse in BAM)
- Quality: double reverse cancels

**Net effect for R1 mapped reverse by BWA (flag=16):**
- PBAT flips strand: `true -> false` (now forward)
- No second revcomp (not negative strand anymore)
- SEQ = revcomp(original)
- Strand flag: `true -> false` (marked forward in BAM)

**e) Tag handling** (line 665-681):
- SAM tags (NM, AS, etc.) are preserved from BWA
- SA and XA tags (supplementary/alternative alignments) are **dropped**

**f) Set pair flags** (line 683-688):
```java
record.setFirstOfPairFlag(!SecondEnd);   // true for R1
record.setReadPairedFlag(true);
record.setSecondOfPairFlag(SecondEnd);   // false for R1
```

### Step 2d: Pairing with `-nonDirectional`

(Bhmem.java line 527-592)

First, reads failing the flag filter (unmapped, not-primary, vendor-fail) are removed. Then, for `noDirectional=true`, **all 4 cross-pair combinations** are tried:

```
Pair A: samsEnd1CT x samsEnd2CT  ->  R1 on GA genome + R2 on GA genome
Pair B: samsEnd1GA x samsEnd2GA  ->  R1 on CT genome + R2 on CT genome
Pair C: samsEnd1CT x samsEnd2GA  ->  R1 on GA genome + R2 on CT genome   (cross-pair)
Pair D: samsEnd1GA x samsEnd2CT  ->  R1 on CT genome + R2 on GA genome   (cross-pair)
```

Without `-nonDirectional` (directional mode), only Pairs A and B would be tried.

For each combination, all R1 hits from one trial are crossed with all R2 hits from the matching trial. Each (R1, R2) candidate must have **opposite strand orientation** (one forward, one reverse).

### Step 2e: Best pair selection

(`comparingSamRecordPbat` method)

Candidate pairs are compared in strict priority order:

| Priority | Criterion | Wins |
|----------|-----------|------|
| 1 | Both mates MAPQ > 0 | Beats either MAPQ = 0 |
| 2 | Sum MAPQ (R1 + R2) | Higher wins |
| 3 | Same chromosome | Same chr beats different chr |
| 4 | Enzyme proximity | Near DpnII site (+-50bp) beats not near |
| 5 | Sum AS (alignment score) | Higher wins |
| 6 | Sum NM (edit distance) | **Lower** wins |
| 7 | Sum CIGAR M-length | Higher wins (longer aligned span) |

**Key**: bhmem **never rejects** for ambiguity. It always picks a winner. If everything is tied, the first candidate encountered is kept.

### Step 2f: Final BAM output

```java
SamPairUtil.setProperPairAndMateInfo(r1, r2, samFileHeader, ORIENTATION, true);
```

Sets mate coordinates, insert size, proper-pair flag. Both FR and RF orientations are accepted. Interchromosomal pairs are output (because `-outputMateDiffChr`).

No MAPQ filter, no dedup -- raw paired BAM output.

---

## Side-by-Side Comparison with Concrete Example

Suppose R1 truly originated from the OB strand and R2 from the OT strand.

### What bismark does:

**R1** (`--pbat`):
```
Original:     ACGTCGATCG
Convert G->A: ACATCAATCA
Trial 1: ACATCAATCA -> CT genome (--nofw, revcomp only) -> CTOT strand
Trial 2: ACATCAATCA -> GA genome (--norc, forward only) -> CTOB strand
-> Pick best by Bowtie2 AS. Reject if tied.
```

**R2** (standard directional):
```
Original:     TGCATGCATG
Convert C->T: TGTATGTATG
Trial 1: TGTATGTATG -> CT genome (--norc, forward only) -> OT strand
Trial 2: TGTATGTATG -> GA genome (--nofw, revcomp only) -> OB strand
-> Pick best by Bowtie2 AS. Reject if tied.
```

R1 and R2 are mapped **independently**. No pairing constraint.

### What bhmem does:

**R1** (pbat):
```
Original:     ACGTCGATCG
Convert G->A: ACATCAATCA
End1CT: ACATCAATCA -> GA genome (pbat swapped)
End1GA: ACATCAATCA -> CT genome (pbat swapped)
```

**R2** (pbat):
```
Original:     TGCATGCATG
Convert C->T: TGTATGTATG
End2CT: TGTATGTATG -> GA genome (pbat swapped)
End2GA: TGTATGTATG -> CT genome (pbat swapped)
```

Then pair all 4 combinations with `-nonDirectional`:
```
Pair A: (R1 -> GA genome) + (R2 -> GA genome)
Pair B: (R1 -> CT genome) + (R2 -> CT genome)
Pair C: (R1 -> GA genome) + (R2 -> CT genome)    <- cross-pair, only with -nonDirectional
Pair D: (R1 -> CT genome) + (R2 -> GA genome)    <- cross-pair, only with -nonDirectional
```

Pick overall best pair by sum MAPQ -> same chr -> enzyme -> AS -> NM -> M-length.

---

## Summary of Equivalences and Differences

### Read conversion: identical

| Mate | Bismark (yap) | Bhmem |
|------|--------------|-------|
| R1 | G->A (`--pbat`) | G->A (`pbat=true`) |
| R2 | C->T (standard directional) | C->T (`pbat=true`) |

Despite different code paths and naming, both tools apply the **same conversion** to each mate.

### Genome trials: identical set

| Mate | Bismark genomes tried | Bhmem genomes tried |
|------|----------------------|---------------------|
| R1 | CT genome + GA genome | GA genome + CT genome |
| R2 | CT genome + GA genome | GA genome + CT genome |

Same two genomes tried for each mate in both tools.

### Key differences

| Aspect | Bismark (yap) | Bhmem |
|--------|--------------|-------|
| **Aligner** | Bowtie2 | BWA mem |
| **R1 mode** | PBAT (`--pbat`) | PBAT (`-pbat`) |
| **R2 mode** | **Standard directional** (no `--pbat`) | **PBAT** (`-pbat` applies to both) |
| **R2 strand restriction** | Bowtie2 `--norc`/`--nofw` restricts orientation per trial | BWA mem reports both orientations freely |
| **Pairing** | Independent single-end per mate | Paired-end: R1+R2 must form opposite-strand pair |
| **Non-directional** | N/A (single-end) | `-nonDirectional` allows cross-genome pairing (Pairs C, D) |
| **Ambiguity** | Rejects reads with tied best AS scores | Never rejects; always picks a winner (first seen if tied) |
| **Best alignment** | Per-read: highest Bowtie2 AS, unique only | Per-pair: sum MAPQ -> same chr -> enzyme -> sum AS -> sum NM -> sum M |
| **Enzyme awareness** | None | DpnII proximity (+-50bp) used in pair selection |
| **Split-read** | Yes (unmapped reads split into fragments, remapped) | No |
| **NM in BAM** | vs **unconverted** genome (includes bisulfite C->T as "mismatches") | vs **converted** reference (BWA tags passed through) |
| **SEQ in BAM** | Original unconverted bases | Original unconverted bases |
| **MAPQ filter** | >= 10 | None (raw MAPQ) |
| **Dedup** | Picard MarkDuplicates | None |
| **Output** | Two independent single-end BAMs merged later | One paired-end BAM |

### Critical difference for R2

Bhmem treats **both R1 and R2 as PBAT**. Yap treats **R1 as PBAT, R2 as standard directional**.

While the read conversions end up the same (R1: G->A, R2: C->T), the **strand/orientation restrictions differ for R2**:

- **Bismark R2** (standard directional): C->T query -> CT genome with `--norc` (forward only, OT), C->T query -> GA genome with `--nofw` (reverse only, OB)
- **Bhmem R2** (PBAT): C->T query -> GA genome and CT genome with **no orientation restriction** from BWA mem

This may explain why R2 shows more discrepancy in cross-pipeline comparisons.

### Effect of `-nonDirectional` in bhmem

Without `-nonDirectional` (directional mode):
- Only Pair A (both from GA genome) and Pair B (both from CT genome)

With `-nonDirectional` (as actually run):
- Also Pair C (R1 from GA, R2 from CT) and Pair D (R1 from CT, R2 from GA)
- This allows R1 and R2 to come from **different** converted genome trials

This is **not equivalent** to bismark's behavior in any mode, since bismark maps each mate independently without cross-genome pairing constraints.
