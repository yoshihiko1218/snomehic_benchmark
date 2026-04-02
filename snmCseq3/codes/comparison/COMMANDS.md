# MapQ / discrepancy benchmark — command log

Paths below use this workspace root:

`/home/jmj7858/epifluidlab/workspace/scnomehic_paper/benchmark/snmCseq3`

Reference FASTA (mm10, must be indexed with `.fai`):

`/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa`

Conda env with `pysam`, `numpy`, `matplotlib`:

`/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python3`

**Maintenance:** When you run new benchmark commands, append a dated subsection under [History](#history) at the bottom of this file (what you ran + one-line purpose).

---

## Environment

```bash
# Optional: activate env
conda activate scnomehic

PY=/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python3
REF=/gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa
cd /home/jmj7858/epifluidlab/workspace/scnomehic_paper/benchmark/snmCseq3
```

---

## 1. Discrepant MAPQ report (yap MAPQ high, bhmem MAPQ low)

Builds per-read TSV with `NM` (bhmem), yap raw/corrected bisulfite mismatch, bhmem CIGAR M/=/X bisulfite counts; summary + multi-page PDF (histograms + NM vs bhmem corrected scatter on page 2).

```bash
$PY codes/comparison/report_discrepant_mismatch.py \
  mapq_comparison/SRR21549292/yap_high_bhmem_low.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  "$REF" \
  -o mapq_comparison/SRR21549292/discrepant_mismatch_report_nmcorr_10k \
  --max-reads 10000
```

- `--max-reads 0` or omit = all rows in the TSV (~679k; long runtime).

---

## 2. Zero–zero location comparison (NM=0 and yap corrected=0)

Compares bhmem vs yap **coordinates** for reads where both mismatch metrics are zero; optional tolerance (bp) and PDF summary.

```bash
$PY codes/comparison/compare_zero_zero_locations.py \
  mapq_comparison/SRR21549292/discrepant_mismatch_report_subset_10k.per_read.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  -o mapq_comparison/SRR21549292/zero_zero_location_compare.tsv \
  --pdf mapq_comparison/SRR21549292/zero_zero_location_report.pdf \
  --tolerance 5
```

- Input must be a `discrepant_mismatch_report*.per_read.tsv` (same columns as produced by §1).

---

## 3. Example read PDFs (SAM-style) for zero–zero categories

Picks up to `-k` read ends per **status** from `zero_zero_location_compare.tsv`, requires `--per-read-tsv` for MAPQ filter + mismatch header; prefers paired fragments by default.

```bash
$PY codes/comparison/make_zero_zero_example_pdf.py \
  mapq_comparison/SRR21549292/zero_zero_location_compare.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  --per-read-tsv mapq_comparison/SRR21549292/discrepant_mismatch_report_subset_10k.per_read.tsv \
  -o mapq_comparison/SRR21549292/zero_zero_example_reads.pdf \
  --per-condition 2
```

---

## 4. Sanity check: aligned-pairs vs CIGAR M/=/X bisulfite counts (bhmem)

```bash
cd codes/comparison
$PY compare_bhmem_bisulfite_methods.py "$REF" ../04.bhmem_bam/SRR21549292.bhmem.bam --max-reads 5000
```

---

## 5. Recompute “NM-like” distance vs `NM:i` (with / without bisulfite mask on subs)

Compares **genomic** edit distance (subs + indels) to the **`NM`** tag; optional TSV and **PDF** scatter (two panels).

```bash
$PY codes/comparison/compare_nm_recompute.py \
  "$REF" \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --max-reads 10000 \
  -o mapq_comparison/SRR21549292/nm_recompute_compare_10k.tsv \
  --pdf mapq_comparison/SRR21549292/nm_recompute_compare_plots.pdf
```

Optional: `--bisulfite-read2-mode pbat_read2` uses strand masking on read 1 and symmetric C/T + G/A masking on read 2 (PBAT / non-directional); default is `strand`.

---

## 6. Standalone bisulfite mismatch CLI (single BAM)

```bash
$PY codes/comparison/bisulfite_corrected_mismatch.py \
  "$REF" \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --method aligned_pairs \
  --max-reads 1000

$PY codes/comparison/bisulfite_corrected_mismatch.py \
  "$REF" \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --method cigar_mx \
  --max-reads 1000
```

- `aligned_pairs`: `get_aligned_pairs(matches_only=True)`.
- `cigar_mx`: explicit **M/=/X** walk (good for heavy soft-clip / indel CIGARs).

---

## 7. NM vs bhmem **corrected substitution** scatter (from per-read TSV)

One-off plot when you already have `discrepant_mismatch_report*.per_read.tsv` with `nm_bhmem` and `bhmem_corrected_mismatch`:

```bash
$PY << 'PY'
import csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
path = "mapq_comparison/SRR21549292/discrepant_mismatch_report_nmcorr_10k.per_read.tsv"
out = "mapq_comparison/SRR21549292/nm_vs_bhmem_corrected_scatter.pdf"
nm, corr = [], []
with open(path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row.get("nm_bhmem", "") == "": continue
        nm.append(float(row["nm_bhmem"])); corr.append(float(row["bhmem_corrected_mismatch"]))
nm, corr = np.array(nm), np.array(corr)
pr = float(np.corrcoef(nm, corr)[0, 1])
fig, ax = plt.subplots(figsize=(6.5, 6))
ax.scatter(nm, corr, alpha=0.22, s=12, rasterized=True)
mx = max(nm.max(), corr.max())
ax.plot([0, mx], [0, mx], "k--", alpha=0.35)
ax.set_xlabel("NM (bhmem)"); ax.set_ylabel("Bhmem corrected mismatch (M/=/X)")
ax.set_title(f"Pearson r = {pr:.4f} (n={len(nm)})")
fig.tight_layout(); fig.savefig(out, dpi=150)
print("Wrote", out)
PY
```

---

## 8. NM recompute **by mate** (read1 vs read2) + `pbat_read2` mask + **MD-based** recompute

3×3 PDF: rows 1–2 = FASTA walks (genomic / strand / `pbat_read2`); row 3 = **MD-based** subs + indels vs `NM:i` (tracks aligner `NM` because `MD` uses the same bisulfite-consistent reference as bhmem, not raw genomic `C` at converted sites). Console prints stats including `from_md`. TSV adds `recompute_from_md`.

```bash
$PY codes/comparison/report_nm_recompute_by_mate.py \
  "$REF" \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --max-reads 12000 \
  -o mapq_comparison/SRR21549292/nm_recompute_by_mate.per_read.tsv \
  --pdf mapq_comparison/SRR21549292/nm_recompute_by_mate_report.pdf
```

Helper in code: `count_nm_style_edit_distance_from_md(read)` in `bisulfite_corrected_mismatch.py` (no FASTA).

---

## 9. Older three-way PDF (bhmem + yap Bowtie2 + yap Bowtie1)

```bash
$PY codes/comparison/make_discrepant_pdf_three.py \
  mapq_comparison/SRR21549292/yap_high_bhmem_low.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  /path/to/yap_bowtie1.bam \
  -o mapq_comparison/SRR21549292/discrepant_reads_three.pdf \
  -n 10
```

---

## 10. **Recommended:** one NM-style recompute for both pipelines (``MD`` + CIGAR)

To **match bhmem ``NM:i`` as closely as possible**, then apply the **identical** computation on yap:

- Use ``recompute_nm_style_from_md`` in ``bisulfite_corrected_mismatch.py`` (alias of ``count_nm_style_edit_distance_from_md``): decode subs from ``MD:Z`` + indel lengths from CIGAR. **No reference FASTA**, no yap-specific vs bhmem-specific branches.
- **Bhmem:** typically very high concordance with ``NM`` (especially read 1); small residuals can appear on read 2 when pooling.
- **Yap/Bismark:** recomputed value usually **equals** the recorded ``NM:i`` on each read.

Cross-pipeline plots / TSV from joined reads:

```bash
$PY codes/comparison/report_nm_md_cross_pipeline.py \
  mapq_comparison/SRR21549292/mapq_comparison.joined.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  -o mapq_comparison/SRR21549292/nm_md_cross_pipeline_50k \
  --max-rows 50000
```

**Contrast:** §11 is a **non-``MD``** genomic walk for a different goal (fair comparison on raw mm10 without using stored ``MD``); it does **not** track bhmem ``NM`` as tightly as this section.

---

## 11. Unified **non-MD** genomic edit distance (same rule on bhmem + yap)

Uses `count_cross_pipeline_comparable_edit_distance` in `bisulfite_corrected_mismatch.py`: **mm10.fa** + CIGAR walk, with **XR-based** masking for Bismark/yap and **strand + read2-symmetric** masking for bhmem when `XR` is absent. **No `MD` tag** — one definition so yap and bhmem are comparable on the **same** reference sequence.

- **Do not** expect yap **`NM:i`** to match this metric (Bismark `NM` is uncorrected vs raw genome).
- Bhmem **`NM:i`** is only **approximated** here (often closer on read 1); for **tag validation** on bhmem, prefer §8 / `count_nm_style_edit_distance_from_md`.

Input: `mapq_comparison.joined.tsv` (`base_id`, `is_r1`, …). Outputs: `*.per_read.tsv`, `*.summary.txt`, `*.plots.pdf`.

```bash
$PY codes/comparison/report_unified_genomic_edit_cross_pipeline.py \
  mapq_comparison/SRR21549292/mapq_comparison.joined.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  "$REF" \
  -o mapq_comparison/SRR21549292/unified_genomic_edit_cross_50k \
  --max-rows 50000
```

---

## History (append new runs below)

### 2026-03-25

- Added this file and `--pdf` to `compare_nm_recompute.py` (NM vs recompute **no bisulfite mask** / **with bisulfite mask**, two-panel scatter + dashed y = x).
- Documented commands used for: `report_discrepant_mismatch`, `compare_zero_zero_locations`, `make_zero_zero_example_pdf`, `compare_bhmem_bisulfite_methods`, `compare_nm_recompute`, `bisulfite_corrected_mismatch`, and inline NM vs bhmem-corrected scatter.

**NM recompute + plot (2k reads, test):**

```bash
/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python3 \
  codes/comparison/compare_nm_recompute.py \
  /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --max-reads 2000 \
  --pdf mapq_comparison/SRR21549292/nm_recompute_compare_plots.pdf
```

Output: `mapq_comparison/SRR21549292/nm_recompute_compare_plots.pdf` (two scatter panels vs `NM:i`).

**NM gap investigation (8k primary reads, residual vs C/T-like subs):**

```bash
/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python3 \
  codes/comparison/investigate_nm_gap.py \
  /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --max-reads 8000 \
  --pdf mapq_comparison/SRR21549292/nm_gap_investigation.pdf
```

Output: `mapq_comparison/SRR21549292/nm_gap_investigation.pdf`; console stats quantify how often naive genomic recompute exceeds `NM` and Pearson(residual, `n_sub_ct`).

**NM recompute by mate + `pbat_read2` in `count_nm_style_edit_distance` (12k primaries):**

```bash
/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python3 \
  codes/comparison/report_nm_recompute_by_mate.py \
  /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  --max-reads 12000 \
  -o mapq_comparison/SRR21549292/nm_recompute_by_mate.per_read.tsv \
  --pdf mapq_comparison/SRR21549292/nm_recompute_by_mate_report.pdf
```

Outputs: `nm_recompute_by_mate_report.pdf`, `nm_recompute_by_mate.per_read.tsv`; added `--bisulfite-read2-mode` to `compare_nm_recompute.py`.

**Cross-pipeline consistent ``NM`` (``MD`` + CIGAR indels) for bhmem + yap:**

Same definition for both BAMs: ``count_nm_style_edit_distance_from_md`` in ``bisulfite_corrected_mismatch.py``. Matches each aligner’s ``NM:i`` closely; use joined MAPQ TSV to pair reads.

```bash
$PY codes/comparison/report_nm_md_cross_pipeline.py \
  mapq_comparison/SRR21549292/mapq_comparison.joined.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  -o mapq_comparison/SRR21549292/nm_md_cross_pipeline_50k \
  --max-rows 50000
```

**Yap corrected mismatch (XR-based) vs bhmem (XR vs mapping strand), 10k discrepant reads:**

```bash
/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python3 \
  codes/comparison/report_discrepant_mismatch.py \
  mapq_comparison/SRR21549292/yap_high_bhmem_low.tsv \
  04.bhmem_bam/SRR21549292.bhmem.bam \
  alignment/Group22/bam/SRR21549292.3C.sorted.bam \
  /gpfs/projects/b1198/epifluidlab/yoshii/reference/mm10/mm10.fa \
  -o mapq_comparison/SRR21549292/discrepant_mismatch_report_nmcorr_10k_xr \
  --max-reads 10000
```

Outputs: `discrepant_mismatch_report_nmcorr_10k_xr.plots.pdf`, `discrepant_mismatch_report_nmcorr_10k_xr.summary.txt`, `discrepant_mismatch_report_nmcorr_10k_xr.per_read.tsv`.

**Cross-pipeline MD-based NM (50k joined rows):** added `report_nm_md_cross_pipeline.py`; ran with `--max-rows 50000` → `nm_md_cross_pipeline_50k.{per_read.tsv,summary.txt,plots.pdf}` (yap: 100% tag==MD-recompute; bhmem pooled R1+R2: ~94.9% exact). This is the **same function** on both BAMs (`recompute_nm_style_from_md`) — the workflow “calibrate to bhmem NM, then apply on yap” without FASTA heuristics.

**Bhmem calibration (8k primaries, local benchmark):** `recompute_nm_style_from_md` vs `NM:i` — pooled ~94.8% exact, mean |Δ| ~0.31; read1 ~99.98% exact; read2 ~89.4% exact. FASTA + `bisulfite_correct=True` (`strand` or `pbat_read2`) stayed far from pooled `NM` (~48% exact, mean |Δ| ~3–4 on the same reads).

**Unified non-MD cross-pipeline edit (50k joined rows):** `report_unified_genomic_edit_cross_pipeline.py` → `unified_genomic_edit_cross_50k.{per_read.tsv,summary.txt,plots.pdf}`. Example summary on this sample: bhmem `NM` vs unified exact match ~52%, mean |Δ| ~3; yap `NM` vs unified ~50% exact, mean |Δ| ~11 (expected). Pearson(unified_bhmem, unified_yap) was ~0 on this slice — interpret together with alignment agreement (different loci/CIGAR for the same read name will decorrelate any per-read metric).
