# Yap vs Bhmem: NM metric, trial picking, MAPQ — handoff for Claude Code

This document summarizes **design decisions, limitations, and implementation pointers** from a working thread on comparing yap and Bhmem on the **same reads** using a **shared NM-style** distance (CIGAR + bisulfite-**converted** reference, CT/GA FASTAs), without unconstrained whole-genome realign as the default goal.

---

## 1. Goal (comparison)

- Compare **yap** vs **bhmem** on the **same reads** with a **shared metric**: recomputed **NM-style** edit distance from **stored CIGAR** + **bisulfite-converted** FASTA (PBAT: CT/GA genome × query-orientation variants), **not** yap’s raw `NM` vs unconverted genomic reference.
- Prefer staying at the **pipeline’s placement** (especially for yap): focus on **path rescoring**, **Bhmem-faithful tie-breaks**, optional **regional `bwa mem`**, rather than full WGS realign.
- Cursor rule (if present): `.cursor/rules/snmCseq3-yap-bhmem-comparison.mdc` scoped to `snmCseq3/codes/comparison/**` and `snmCseq3/mapq_comparison/**`.

---

## 2. What Bhmem actually does (single-end)

Bhmem’s single-record merge (`comparingSamRecord`) orders candidates by:

1. **MAPQ** (higher wins)  
2. **AS** (higher wins)  
3. **NM** (lower wins)  
4. **CIGAR `M` span** (Bhmem/htsjdk: count **only** op `M`, not `=`/`X`)

Paired PBAT uses `comparingSamRecordPbat` (sum MAPQ, enzyme tie-break if loaded, sum AS, sum NM, sum `M`, etc.). See `bhmem_equivalent_selection.py` and comments referencing `Bhmem.java`.

---

## 3. Fixed BAM: why MAPQ does not separate four trials

On **one merged primary alignment**, you have **one** `MAPQ`, **one** `POS`+`CIGAR`, usually **one** `AS:i` from the winning line.

The four **PBAT trials** (CT vs GA × two converted-query shapes) are **counterfactual interpretations** of the **same** placement. They share:

- the same **MAPQ**  
- the same **CIGAR** (hence same `M` span)  
- usually the same **AS** tag  

So applying the **full** Bhmem order is **faithful**, but **MAPQ first is a no-op** for trial choice: it never breaks ties. The first field that typically **differs** across trials is **recomputed NM** (per trial on converted FASTA with fixed CIGAR). **True** per-trial **AS** would differ only if you had **separate BWA lines** per pass (e.g. jbwa / full CT vs GA realign).

**You cannot recover four “original” MAPQs** from the merged BAM: MAPQ was computed for the **winning** hit in the **full** index search; the other three trials’ MAPQs were never stored.

---

## 4. Regional `bwa mem` on a window: why local MAPQ/AS/NM is not “the same” criterion

Regional BWA (slice of converted FASTA around the locus) produces **MAPQ/AS/NM for that rerun** against a **truncated reference**.

- **Local MAPQ** = uniqueness **inside the window**, not genome-wide; the set of competing loci is **wrong** relative to the original run.  
- **AS/NM** refer to **where BWA places the read in that run**; they need not match **`NM:i` / `AS:i`** on the **stored** line if CIGAR/pos differ slightly.  
- **Anchoring** (`anchor_slop` near BAM start) reduces off-locus hits but does not reproduce full Bhmem.

Empirical sketch from prior benchmarks (e.g. `SRR21549289.bhmem.bam`, pad ~6000): regional pick vs `NM:i` was **~87%** vs **~88.6%** for **min distance over four path walks** on fixed CIGAR; **`NM` in one of the four trial distances** was very high (~**99%**). So pain is **trial disambiguation** when multiple trials look good, not “NM nowhere near the set.”

**Bottom line:** window MAPQ can be a **heuristic** tie-breaker; it is **not** a substitute for **original** MAPQ for picking the same alignment Bhmem would have reported.

---

## 5. MD vs `NM:i` on Bhmem BAM

**MD-based** recompute often tracks **`NM:i`**, but Bhmem can **rebuild `SEQ`** from FASTQ while **`NM`/`MD` stay from BWA’s SAM line** → occasional **SEQ vs tag** inconsistency (~few percent in samples). So MD recompute ≠ `NM:i` in every read even when both are well-defined.

---

## 6. Implementation map (repository)

| Piece | Role |
|--------|------|
| `bhmem_equivalent_selection.py` | Bhmem-style single/pair ordering; `EnzymeRegionIndex`; `enumerate_pbat_single_trials`; `recompute_nm_bhmem_style_single_pbat`; **`pick_pbat_single_trial_margin_at_min_dist`**; **`pick_pbat_single_trial`** (`bhmem_fold`, `margin_at_min_dist`, NM-guided strategies). |
| `bisulfite_corrected_mismatch.py` | `count_nm_style_edit_distance_converted_explicit`, `pbat_converted_genome_trial_distances`, `recompute_nm_from_converted_genomes_pbat_no_md` (delegates to Bhmem-style single fold; on fixed BAM collapses to **min trial dist** when MQ/AS/M tied). |
| `regional_bwa_trial_pick.py` | Per-trial regional `bwa index` + `bwa mem`, parse primary MAPQ/AS/NM/M, fold with `bhmem_prefer_second_single`; `--anchor-slop`. |
| `validate_bhmem_nm_recompute.py` | MD vs `NM:i`, four-trial stats, pair fold vs tags. |
| `benchmark_trial_pick_strategies.py` | Compare pick strategies vs `NM:i` / unique-NM-trial label. |
| `cigar_nm_walk/` | Fast fixed-CIGAR NM walk (C). |

**Trial disambiguation strategies (tag-free vs calibration):**

- **`margin_at_min_dist`:** among global **min** recomputed distance, maximize **margin** `min(other dist) - dist`; then Bhmem fold on ties. No `NM:i` required.  
- **`unique_nm_match_else_*` / `nearest_nm`:** use **`NM:i`** to disambiguate — useful **only** when the tag is the reference (Bhmem calibration); **circular** if defining “truth” as the same tag.

---

## 7. Answers to explicit questions (short)

| Question | Answer |
|----------|--------|
| Can you get original MAPQ for all four trials from the BAM? | **No** — only the winner’s MAPQ is stored; others need full per-pass realign or saved candidates. |
| Does window BWA MAPQ correctly pick the “original” alignment? | **Not as a strict criterion** — it’s MAPQ in a **different** search space; heuristic at best. |
| Should you still use MAPQ first like the original code? | **Yes** for **faithful comparator** between full candidates; on **fixed** BAM four-way scoring it’s **correct but usually tied**, so NM (recomputed) or per-trial AS does the real separation. |
| What to use to “retrieve NM” fairly for cross-pipeline comparison? | **Recomputed NM-style** on **converted** FASTA with **fixed** CIGAR, with trial choice via full tuple where fields differ; **`NM:i` on yap** is often not the right target for bisulfite-aware comparison. |

---

## 8. Suggested next steps (if extending)

- **Mate-consistent** joint assignment (both reads share conversion assignment) where single-end is ambiguous.  
- **jbwa / per-trial full align** to obtain real **per-trial AS** (and MAPQ) in the same search space as Bhmem.  
- Document in paper: **fixed-CIGAR min/margin** vs **tag `NM:i`** agreement limits and **regional BWA** as approximate only.

---

*Generated as a conversation handoff; adjust paths and numbers if your environment or BAM differs.*
