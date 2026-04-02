# Bug in `Bhmem.java` `String2SamRecord`: `originalSeq` mutation across SAM lines

## Location

`/home/jmj7858/epifluidlab/software/bisulfitehic/src/main/java/edu/mit/compbio/bisulfitehic/mapping/Bhmem.java`

Method: `String2SamRecord` (line 620)

## The problem

When BWA outputs **multiple SAM lines** for a single read (primary + supplementary for chimeric/split reads), the method loops over them:

```java
for(String line : lines.split("\\n")){
    // ... uses and mutates `originalSeq` ...
}
```

The variable `originalSeq` is a method parameter (passed by reference as a Java `String`, but reassigned within the loop body). The PBAT strand-handling block **reassigns** `originalSeq` in each iteration:

```java
// Line 648-653
if(pbatThisEnd){
    record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
    originalSeq = SequenceUtil.reverseComplement(originalSeq);  // <-- MUTATES
    baseQ = new StringBuilder(baseQ).reverse().toString();
}

// Line 656-660
if(record.getReadNegativeStrandFlag()){
    originalSeq = SequenceUtil.reverseComplement(originalSeq);  // <-- MUTATES AGAIN
    baseQ = new StringBuilder(baseQ).reverse().toString();
}

record.setReadString(modifySeqByCigar(originalSeq, record.getCigarString()));
```

After iteration 1 finishes, `originalSeq` has been reverse-complemented 0, 1, or 2 times depending on BWA's strand flag for that line. Iteration 2 starts with this **already-modified** `originalSeq`, not the original FASTQ sequence. The PBAT flip and strand check are then applied again on wrong input.

## Step-by-step corruption trace

Let `orig` = original FASTQ sequence, `RC(x)` = reverse complement of x.

### Case: Both BWA lines map to reverse strand

**Iteration 1** (BWA flag: reverse strand):
1. PBAT block: flip flag (now forward), `originalSeq = RC(orig)`
2. Strand check: flag is now forward → skip
3. SEQ stored = `modifySeqByCigar(RC(orig), cigar1)` → **correct** (RC for reverse-mapped PBAT)
4. **`originalSeq` is now `RC(orig)` entering iteration 2**

**Iteration 2** (BWA flag: reverse strand):
1. PBAT block: flip flag (now forward), `originalSeq = RC(RC(orig)) = orig`
2. Strand check: flag is now forward → skip
3. SEQ stored = `modifySeqByCigar(orig, cigar2)` → **WRONG**, should be `RC(orig)`

### Case: Iteration 1 reverse, iteration 2 forward

**Iteration 1** (reverse):
1. PBAT: flip→forward, `originalSeq = RC(orig)`
2. Strand: forward → skip. `originalSeq` = `RC(orig)`

**Iteration 2** (forward):
1. PBAT: flip→reverse, `originalSeq = RC(RC(orig)) = orig`
2. Strand: reverse → `originalSeq = RC(orig)`
3. SEQ = `modifySeqByCigar(RC(orig), cigar2)` → **WRONG**, should be `orig` (forward-mapped PBAT)

### Cases that happen to work

- Both forward: each iteration does PBAT RC then strand RC, net = identity. `originalSeq` returns to `orig` between iterations. (**Accidentally correct.**)
- Iteration 1 forward, iteration 2 reverse: similar accidental correctness.

## Which reads are affected

- Only reads where BWA outputs **2+ SAM lines** (chimeric/split/supplementary alignments)
- The corruption depends on the **strand combination** between lines
- These reads typically have **heavy soft-clipping** (e.g. `76S27M28S`) and **low MAPQ** (0 or single digits)
- Affects ~5% of R2 reads, ~0.7% of R1 reads in tested samples

## Observable effect

- The BAM stores **wrong SEQ** for the affected alignment
- **NM:i** and **MD:Z** tags are **correct** (they come from BWA's original SAM output, computed before Bhmem rebuilds SEQ)
- When you try to recompute NM by walking the stored SEQ against the converted reference, you get a completely different value (e.g., NM tag says 0 but recomputed trials give [18, 20, 20, 20])
- Downstream bisulfite calling (allcools, Bis-SNP) reads the **stored SEQ** and may call incorrect methylation at those positions

## Suggested fix

Save a local copy of `originalSeq` at the start of each loop iteration:

```java
for(String line : lines.split("\\n")){
    String iterSeq = originalSeq;   // fresh copy each iteration
    String iterBaseQ = baseQ;       // fresh copy of base quality too

    // ... all subsequent code uses iterSeq instead of originalSeq ...

    if(pbatThisEnd){
        record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
        iterSeq = SequenceUtil.reverseComplement(iterSeq);
        iterBaseQ = new StringBuilder(iterBaseQ).reverse().toString();
    }
    if(record.getReadNegativeStrandFlag()){
        iterSeq = SequenceUtil.reverseComplement(iterSeq);
        iterBaseQ = new StringBuilder(iterBaseQ).reverse().toString();
    }
    record.setReadString(modifySeqByCigar(iterSeq, record.getCigarString()));
    record.setBaseQualityString(iterBaseQ);
    // ...
}
```

Note: `baseQ` has the exact same mutation problem and should also be copied per iteration.

## Scope of impact

- **NM/MD tags**: Unaffected (sourced from BWA SAM output)
- **SEQ in BAM**: Wrong for ~5% of R2 reads (chimeric/split reads with certain strand combinations)
- **Base qualities**: Also potentially wrong (same mutation pattern on `baseQ`)
- **Methylation calling**: Any caller reading the stored SEQ for affected reads will see wrong bases → wrong methylation calls
- **Recomputation from BAM**: Cannot recover NM from CIGAR + converted FASTA for affected reads; must trust the NM tag

## Validation

Tested on 3 samples (SRR21549289, SRR21549292, SRR21549298), 3000-5000 reads each. Reads where NM:i is not found in any of the 4 converted-genome trial distances are consistent with the corruption pattern: almost all have heavy soft-clipping, low MAPQ, NM=0, and trial distances >>0.
