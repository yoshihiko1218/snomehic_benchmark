# snmCAT / snmC2T-seq (NOMe) pipeline — final runbook

End-to-end yap `mct --nome` pipeline, validated on brain snmC2T-seq (UMB5580):
GCH ~16–19% >> HCH ~4–5% = real chromatin accessibility. All bug fixes are baked into the scripts.

Conda env: `mapping` (has yap 1.6.9, bismark, STAR 2.7.3a, allcools).
BASE = `/gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT`

## One-time prerequisites (already done — reuse)
- STAR index matching the env's STAR (2.7.3a):
  `/gpfs/projects/b1198/.../reference/hg38/star_2.7.10a_gencode.v36_sjdb100`  (built by `codes/04.build_star_index_2.7.10a.sh`)
- NOMe mapping config `codes/mapping_config_nome.ini`, generated with:
  ```
  yap default-mapping-config --mode mct --barcode_version V2 --nome \
    --bismark_ref   <ref>/hg38_bismark \
    --genome_fasta  <ref>/hg38_bismark/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --star_ref      <ref>/hg38/star_2.7.10a_gencode.v36_sjdb100 \
    --gtf           <ref>/hg38/gencode.v36.annotation.gtf \
    --chrom_size_path <ref>/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes
  ```

## Per-batch pipeline (brain snmC2T-seq = codes 15–19)

```bash
cd /gpfs/projects/b1042/epifluidlab/yoshii/scnomehic_paper/benchmark/snmCAT
source ~/.bashrc && conda activate mapping

# 0. Build the download list: pick a snmC2T-seq batch from meta.tsv (col4 = "R1ftp;R2ftp")
#    (already done -> codes/download_list_brain.txt: 100 cells, 190321_mCTseq_hs_29yr / UMB5580)
awk -F'\t' 'NR>1 && index($3,"190321_mCTseq_hs_29yr"){print $4; n++} n>=100{exit}' meta.tsv > codes/download_list_brain.txt

# 1. Download                -> fastq_brain/
sbatch codes/15.download_brain.sh                 # SLURM array 1-100%10
#    wait until: ls fastq_brain/*.fastq.gz | wc -l  == 200

# 2. yap-pattern symlinks
bash codes/16.rename_symlink_brain.sh             # -> *-R[12].fq.gz

# 3. Generate per-cell snakemake (NOMe config)
yap start-from-cell-fastq -o mapping_brain \
    -config codes/mapping_config_nome.ini \
    -fq "fastq_brain/*-R[12].fq.gz"

# 4. Patch Snakefiles (yap leaves these blank — must inject):
#    bismark_reference, star_reference, nome_flag_str='--nome'
bash codes/17.patch_nome_brain.sh

# 5. Map (yap mct --nome): one array task per Group
N=$(grep -cve '^\s*$' mapping_brain/snakemake/snakemake_cmd.txt)
sbatch --array=1-$N codes/18.run_brain_nome.sh    # bismark + STAR + select-dna --nome + allc
#    wait until: ls mapping_brain/Group*/allc/*.allc.tsv.gz | wc -l == 100

# 6. QC summary             -> mapping_brain/stats/MappingSummary.csv.gz
yap summary -o mapping_brain

# 7. Per-cell HCG/GCH/HCH counts + methylation rates
sbatch codes/19.collect_brain_gch.sh              # -> mapping_brain/stats/hcg_gch_nome.tsv
```

## Outputs
- `mapping_brain/Group*/allc/<cell>.allc.tsv.gz` — per-cytosine methylation, NOMe context (upstream+C+2down).
- `mapping_brain/Group*/rna_bam/*.feature_count.tsv` — RNA gene counts.
- `mapping_brain/stats/MappingSummary.csv.gz` — 92-col per-cell QC.
- `mapping_brain/stats/hcg_gch_nome.tsv` — cell_id, HCG/GCH/HCH site counts + mC rates.
  - HCG (=H-CG) = CpG methylation. GCH (=G-CH) = NOMe accessibility. GCG excluded.

## Key fixes already in the scripts (why a naive yap run fails)
1. STAR index must match the env STAR version (2.7.3a) — the 2.7.11b index gives `genomeType` FATAL.
2. `source ~/.bashrc` must come BEFORE `set -u` (bashrc has an unbound var → 0-sec failure).
3. Runner does `--rerun-incomplete --unlock` then `--rerun-incomplete` (stale-lock + IncompleteFiles safe).
4. yap omits `bismark_reference`/`star_reference` and leaves `nome_flag_str` blank → the patch script injects them.
5. NOMe needs both `num_upstr_bases=1` (config) AND `nome_flag_str='--nome'` (select-dna-reads excludes GpC).

## To run a DIFFERENT batch
Swap the title filter in step 0 (e.g. `190305_mCTseq_hs_21yr` for UMB5577), point the dirs/scripts at a new
name, and repeat. Only snmC2T-seq (NOMe) batches show GCH accessibility — see memory `gse140493-nome-batches`.
The H1/HEK `mCT` batches (171009 scmCT, 180615 snmCT) are non-NOMe (GCH ≈ background).
```
