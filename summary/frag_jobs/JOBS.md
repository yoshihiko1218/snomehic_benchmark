# Uniform fragment-level alignment QC — jobs

Goal: recompute metric 2 (uniquely-mapped count) + metric 3 (MapQ30 rate) for ALL
datasets with ONE consistent definition (decided 2026-06-30):
- **unit = fragment** (per molecule = unique read-name; R1+R2 collapse, no double-count)
- **duplicates excluded** everywhere
- **rate = MapQ30 fragments / uniquely-mapped fragments** (MapQ30 / mapped)
- a fragment is uniquely-mapped if it has >=1 primary, mapped, non-duplicate read (OR);
  MapQ30 if >=1 such read has MAPQ>=30.

Code: `summary/frag_counts.py` (per-cell), `summary/frag_jobs/frag_array.sh` (array,
markdup-if-needed), `summary/frag_droplethic.py` + `frag_droplethic.sbatch` (per-barcode).
Manifests: `summary/frag_jobs/<ds>.manifest.tsv` (built by make_manifests.py).
Outputs: `summary/frag_counts/<ds>/<cell>.tsv` (cols cell, uniqmap_frag, mapq30_frag, rate).

## BAM source per dataset (the "correct BAM")
| ds | BAM | dedup | markdup step? |
|----|-----|-------|----|
| nagano | alignment/{cell}.markdup.bam | flagged | no |
| smallwood | 05.align_mm10/{cell}.rmdup.RG.bam | flagged | no |
| scnome | 04.alignment/{cell}_{1,2}.rmdup.RG.bam (R1+R2) | flagged | no |
| snmCseq3 | 04.bhmem_bam/{cell}.calmd.bam | flagged | no |
| scnomehic | EXTERNAL gm_sc_new/{cell}.calmd.bam | flagged | no |
| snmCseq2 | **yap_mapping_mm10/Group*/bam/{cell}.final.bam** (YAP, deduped) | removed | no |
| snmCAT | mapping_brain/Group*/bam/{cell}.dna_reads.bam | none | **YES** |
| droplethic | my_project/03.mapping/SRR27586278_hg38.bam (CB-tagged, per-barcode) | flagged | no |

## Submitted jobs (2026-06-30, account b1042 / genomics)
First batch failed (set -u before sourcing ~/.bashrc -> BASHRCSOURCED unbound). Fixed
(moved `set -uo pipefail` after conda activate). Resubmitted:
| ds | jobid | array |
|----|-------|-------|
| nagano | 5530301 | 1-15 |
| smallwood | 5530302 | 1-51 |
| scnome | 5530303 | 1-23 |
| snmCseq3 | 5530304 | 1-100 |
| scnomehic | 5530305 | 1-187 |
| snmCseq2 | 5530306 | 1-96 (YAP final.bam) |
| snmCAT | 5530307 | 1-100 (markdup first) |
| droplethic | 5530308 | single (198G scan, ~hours) |

Check: `squeue -u jmj7858 | grep frag`
Collect when done: `python summary/frag_jobs/collect.py` -> summary/frag_counts_all.tsv

## UPDATE 2026-06-30 (later): compute BOTH before-dedup and after-dedup
User: calculate both before-dedup (include dups) and after-dedup (exclude dups) for
count + rate. frag_counts.py / frag_droplethic.py now emit 6 values per cell:
uniq_preDedup, mapq30_preDedup, rate_preDedup, uniq_postDedup, mapq30_postDedup, rate_postDedup.
snmCseq2 REVERTED from YAP final.bam (dups already removed -> no before-dedup) to raw
clean_bismark_bt2.bam + markdup (dups flagged -> both versions). Resubmitted:
nagano 5531418, smallwood 5531419, scnome 5531420, snmCseq3 5531421, scnomehic 5531422,
snmCseq2 5531423, snmCAT 5531424, droplethic 5531425.
Example nagano SRR921526: rate_preDedup 85.5% vs rate_postDedup 63.9% (dups are MapQ30-rich).

## UPDATE 2026-06-30: snmCseq3 contacts at BOTH MapQ10 and MapQ30
Correction: YAP CisLong = cis>1kb (min_gap=1000), NOT 2.5kb (earlier claim was wrong).
YAP snmCseq3 contacts are deduped + 1kb, but filtered at MapQ>=10 (Snakefile
`samtools view -q 10`), while the other Hi-C methods use MapQ>=30. So keep BOTH:
- MapQ10 (YAP native): from alignment/stats/MappingSummary.csv.gz
  (CisShortContact/CisLongContact/TransContact) -- already available.
- MapQ30: re-filter YAP 3C.sorted.bam at q30 and re-run `yap-internal generate-contacts`
  (mapping env; --min_gap 1000). Job 5541098 (array 1-100) ->
  summary/frag_counts/snmCseq3_contacts_q30/<cell>.q30.counts.txt (CisShort,CisLong,Trans).
One-cell validation (SRR21549383): q10 cis>1kb=8946 trans=10056 ; q30 cis>1kb=7177 trans=7823.
cis_n = CisShort+CisLong ; cis_gt1kb = CisLong ; trans_ratio = Trans/(CisShort+CisLong+Trans).

## NEXT (after jobs finish)
Replace the metric-2 (readcount) and metric-3 (mapq30) sources in qc.ipynb with these
uniform fragment counts, re-run notebook. Conversion figures + loci unaffected.
