# snmCseq3 MAPQ Comparison: yap (Bowtie2) vs bhmem (BwaMem)

Compare mapping quality (MAPQ) of the **same reads** after alignment with two different pipelines:

- **yap**: Bismark + Bowtie2
- **bhmem**: BwaMem-based bisulfite aligner

## Files Used

| Pipeline | BAM location | Notes |
|----------|--------------|-------|
| bhmem | `04.bhmem_bam/{sample}.bhmem.bam` | Raw BwaMem output; MAPQ from aligner |
| yap | `alignment/Group{N}/bam/{sample}.3C.sorted.bam` | Merge of R1+R2 two_mapping (Bowtie2); MAPQ preserved |

Both pipelines start from the same trimmed FASTQ (cutadapt with same parameters). Read IDs are matched by `(base_id, is_r1)`:
- bhmem: QNAME = `SRR21549292.N`, FLAG distinguishes R1/R2
- yap: QNAME = `SRR21549292.N_1_1` or `SRR21549292.N_2_2` (N_1=R1, N_2=R2)

## Usage

### Single sample (default: SRR21549292)

```bash
./run_mapq_comparison.sh
```

Or for a specific sample:

```bash
./run_mapq_comparison.sh SRR21549298
```

### Manual steps

```bash
# 1. Extract MAPQ from each BAM
python extract_mapq.py bhmem 04.bhmem_bam/SRR21549292.bhmem.bam -o bhmem_mapq.tsv
python extract_mapq.py yap   alignment/Group22/bam/SRR21549292.3C.sorted.bam -o yap_mapq.tsv --primary-only

# 2. Compare and optionally plot
python compare_mapq.py bhmem_mapq.tsv yap_mapq.tsv -o mapq_comparison --plot
```

### Run on all overlapping samples

```bash
./run_mapq_comparison_batch.sh
```

## Output

- `{sample}/bhmem_mapq.tsv`, `yap_mapq.tsv`: Extracted (base_id, is_r1, mapq)
- `{sample}/mapq_comparison.joined.tsv`: Matched reads with both MAPQs and difference
- `{sample}/mapq_comparison.stats.txt`: Summary (mean, median, correlation, etc.)
- `{sample}/mapq_comparison.plots.pdf`: Scatter, diff histogram, distributions (if matplotlib available)

## Three-way PDF (bhmem + yap Bowtie2 + yap Bowtie1)

`make_discrepant_pdf_three.py` builds one page per read with all three pipelines. Bowtie1 uses `alignment_bowtie1/Group25/` (same sample `SRR21549292` as Bowtie2 `alignment/Group22/`).

```bash
conda activate scnomehic
python codes/comparison/make_discrepant_pdf_three.py \
  mapq_comparison/SRR21549292/yap_high_bhmem_low.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  alignment_bowtie1/Group25/bam/SRR21549292.3C.sorted.bam \
  -o mapq_comparison/SRR21549292/discrepant_reads_three.pdf -n 10
```

Or: `bash codes/comparison/run_discrepant_pdf_three.sh` (expects BAMs under `snmCseq3/`).

## Extract subset BAMs (bhmem + yap Bowtie2 + yap Bowtie1)

Writes **sorted, indexed** BAMs for the first N read `base_id`s from `yap_high_bhmem_low.tsv` (default N=10), same IDs as the PDFs.

```bash
cd snmCseq3
bash codes/comparison/extract_discrepant_subset_bams.sh 10
```

Output: `mapq_comparison/SRR21549292/subset_bams/`

- `discrepant_subset.bhmem.sorted.bam` (+ `.bai`)
- `discrepant_subset.yap_bowtie2.sorted.bam` (+ `.bai`)
- `discrepant_subset.yap_bowtie1.sorted.bam` (+ `.bai`) — skipped if Bowtie1 BAM missing under `alignment_bowtie1/Group25/`
- `discrepant_subset_read_ids.txt`

BAMs are **coordinate-sorted** so `samtools index` works. Requires `samtools` on `PATH`.

### Manual `samtools` commands (same read list for all three)

Create the read-name list once (first 10 `base_id`s from the discrepant set, or edit as needed):

```bash
cat > /tmp/srr21549292_pdf_reads.txt << 'EOF'
SRR21549292.1
SRR21549292.10
SRR21549292.100000
SRR21549292.1000004
SRR21549292.1000005
SRR21549292.1000007
SRR21549292.1000009
SRR21549292.1000011
SRR21549292.1000013
SRR21549292.1000016
EOF
```

**bhmem** — QNAME equals `base_id` exactly:

```bash
cd snmCseq3
samtools view -h -N /tmp/srr21549292_pdf_reads.txt 04.bhmem_bam/SRR21549292.bhmem.bam \
  | samtools sort -o pdf_subset.bhmem.sorted.bam -
samtools index pdf_subset.bhmem.sorted.bam
```

**yap Bowtie2** (`alignment/Group22/`) — QNAME looks like `SRR21549292.1_1_1`, `SRR21549292.1_2_2`, splits `..._1-l`, etc. Filter with `base_id` + `_`:

```bash
samtools view -h alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  | awk -v f=/tmp/srr21549292_pdf_reads.txt '
    BEGIN { while ((getline l < f) > 0) n[l]=1 }
    /^@/ { print; next }
    { for (k in n) if (index($1, k "_") == 1) { print; next } }
  ' \
  | samtools sort -o pdf_subset.yap_bowtie2.sorted.bam -
samtools index pdf_subset.yap_bowtie2.sorted.bam
```

**yap Bowtie1** (`alignment_bowtie1/Group25/`) — **same `awk` filter** as Bowtie2. QNAMEs differ only in naming (`SRR21549292.1_1` / `_2` instead of `_1_1` / `_2_2`); prefix still matches `base_id` + `_`.

```bash
samtools view -h alignment_bowtie1/Group25/bam/SRR21549292.3C.sorted.bam \
  | awk -v f=/tmp/srr21549292_pdf_reads.txt '
    BEGIN { while ((getline l < f) > 0) n[l]=1 }
    /^@/ { print; next }
    { for (k in n) if (index($1, k "_") == 1) { print; next } }
  ' \
  | samtools sort -o pdf_subset.yap_bowtie1.sorted.bam -
samtools index pdf_subset.yap_bowtie1.sorted.bam
```

**Single read pair** example (`SRR21549292.1`):

```bash
# bhmem
samtools view -h 04.bhmem_bam/SRR21549292.bhmem.bam \
  | awk '$1 ~ /^@/ || $1=="SRR21549292.1"' \
  | samtools sort -o one.bhmem.sorted.bam -
samtools index one.bhmem.sorted.bam

# yap Bowtie2
samtools view -h alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  | awk '$1 ~ /^@/ || $1 ~ /^SRR21549292\.1_/' \
  | samtools sort -o one.yap_bowtie2.sorted.bam -
samtools index one.yap_bowtie2.sorted.bam

# yap Bowtie1 (same regex: lines starting with SRR21549292.1_)
samtools view -h alignment_bowtie1/Group25/bam/SRR21549292.3C.sorted.bam \
  | awk '$1 ~ /^@/ || $1 ~ /^SRR21549292\.1_/' \
  | samtools sort -o one.yap_bowtie1.sorted.bam -
samtools index one.yap_bowtie1.sorted.bam
```

Always **`samtools sort` before `samtools index`** on subsets, or indexing will fail.

## Discrepant-read summary (yap>30, bhmem<30)

For reads where yap MAPQ > 30 and bhmem MAPQ < 30, `summarize_discrepant_reads.py` produces:

- **discrepant_summary.summary.txt**: Same chr?, same pos?, position distance (when same chr), mismatches (NM), strand agreement
- **discrepant_summary.detail.tsv**: Per-read chr, pos, NM, distance, etc.

**NM note**: bhmem aligns to bisulfite-converted reference (C→T at CpG are matches), so NM is low. yap/Bismark uses unconverted reference; C-to-T conversions count as mismatches, so NM is higher (~20–40). The two NM values are not directly comparable.

## Why NM differs (same locus: bhmem low, yap/Bismark high)

- **Bismark (yap)** writes SAM against the **genomic** reference. `NM` / `MD` count differences to the **unconverted** genome. After bisulfite treatment, unmethylated **C** in DNA becomes **T** in the read, but the reference still has **C** → those look like “mismatches” and dominate `NM` (e.g. `NM:i:31` vs `NM:i:1` for the same read).
- **Bhmem** aligns to **bisulfite-converted** indexes (BWA-style); expected **C→T** matches the index, so **`NM` stays small**.

Bismark adds **`XM`**, **`XR`**, **`XG`** for methylation context; raw `NM` is still not bisulfite-normalized.

### Bisulfite-corrected mismatch (heuristic)

`bisulfite_corrected_mismatch.py` walks aligned bases (read + `reference.fa`) and counts mismatches **excluding**:

- Forward-mapped reads (`FLAG` without reverse): **C** (ref) vs **T** (read)
- Reverse-mapped reads: **G** (ref) vs **A** (read), complement of C→T on the other strand

This assumes full conversion at unmethylated C; **methylated C** (still C in read after failed conversion or biological mC) is not distinguished here — use methylation tools / `XM` for that.

```bash
conda activate scnomehic
python codes/comparison/bisulfite_corrected_mismatch.py \
  /path/to/mm10.fa \
  mapq_comparison/SRR21549292/subset_bams/discrepant_subset.yap_bowtie2.sorted.bam \
  --max-reads 100
```

## Dependencies

- `samtools` (required)
- `python3` with `numpy` (for compare_mapq)
- `matplotlib` (optional, for plots)
- `pysam` (optional for streaming; **required** for `bisulfite_corrected_mismatch.py`)
