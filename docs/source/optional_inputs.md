# Optional input files

Besides the fastqs and the references, some stages read tables or files that you provide. This page lists each one: which config key points to it, which stage reads it and when, and its format. All paths go in the config; `null` turns an optional input off.

## Overview

| Input | Config key | Read by | Needed when |
|---|---|---|---|
| [Fastq list](#fastq-list) | `experiment: inputs: select_fastqs` | every run from the fastqs | `workflow: start: align` (not read with `start: split_bams` or `wasp`) |
| [Wet-lab sheet](#wet-lab-sheet) | `experiment: inputs: wet_lab_info` | `kite` (HTO counts), `hashsolo` demultiplexing | HTO runs (`experiment: hto: enabled`, `demultiplexing: method: solo` or `both`) |
| [Pool info (vcf_info)](#pool-info-vcf_info) | `modules: gt_demux: vcf_info` | cellSNP, vireo, STARsolo's WASP mode | genotype demultiplexing (`method: vireo` or `both`); `STARsolo: wasp` |
| [All-donor VCF](#all-donor-vcf) | `modules: gt_demux: donors_per_pool_file` | vireo | genotyped runs with a 2-column `vcf_info` |
| [Reference SNP sites](#reference-snp-sites) | `modules: gt_demux: genome_1k_ref` | cellSNP | pileups without a per-pool VCF (e.g. `run_tag: wo_gt`) |
| [Vireo donor-name converter](#vireo-donor-name-converter) | `modules: gt_demux: donorName_conv` | demultiplexing (vireo h5ads) | vireo's donor names need renaming in the h5ad |
| [Swap correction file](#swap-correction-file) | `modules: demultiplex: swap_correction` | final h5ads (`demux_final`), `split_bams` | confirmed donors per pool; `split_by: input: final`, `start: split_bams` |
| [Demultiplex versions](#demultiplex-versions) | `modules: demultiplex: swap_correction: versions` | as the swap correction file | the version column holds names, not `<tool>:<run_tag>` |
| [External donor labels](#external-donor-labels) | `modules: split_bams: split_by: external: file` | `split_bams` | `split_by: input: external` |
| [Pools table (bams.tsv)](#pools-table-bamstsv) | `modules: split_bams: split_by: external: pools` | `split_bams` | `start: split_bams` with bams outside the aligner's folders |
| [Split-bam donor-name converter](#split-bam-donor-name-converter) | `modules: split_bams: donor_name_converter` | `split_bams` | per-donor bams should carry other names |
| [Reference genotypes for MBV](#reference-genotypes-for-mbv) | `modules: identify_swaps: ref_vcf` | `identify_swaps` (QTLtools mbv) | `workflow: until: identify_swaps` |
| [WASP samples table (wasp_samples.tsv)](#wasp-samples-table-wasp_samplestsv) | `modules: wasp: samples` | `wasp` | `workflow: until: wasp` or `start: wasp` |
| [WASP sample names](#wasp-sample-names) | `modules: wasp: sample_names` | `wasp` | donors were swap-corrected (`final_donor`), or pools or donors aren't named as in the genotypes, and the samples table leaves `final_donor` / `vcfID` empty |
| [WASP genotype sets](#wasp-genotype-sets) | `modules: wasp: genotypes` | `wasp` | `workflow: until: wasp` or `start: wasp` |
| [WASP remap references](#wasp-remap-references) | `modules: wasp: rna: remap`, `atac: remap` | `wasp` | `start: wasp`; or a genome other than the aligner's |

General rules for the tables:

- **Separator by extension.** `.tsv` is tab-separated and `.csv` comma-separated. Where `.txt` is accepted, it is noted for that file (the pool info file reads `.txt` as space-separated; the external labels read it as tab-separated).
- **Header row.** Every table has one. Where the columns are found by position instead of name, the section says so.
- **Pool and donor names** must match the `{pool}` names of the run (the fastq list, or the pools of the labels) and can't contain `/`.
- **Paths** in tables are best given as absolute paths: jobs may run from another directory on the cluster.

## Fastq list

`experiment: inputs: select_fastqs`. Required with `workflow: start: align`; not read when the run starts at `split_bams` or `wasp`.

A text file, one pool per line, as its fastq folder name (`-cDNA` is removed to get the `{pool}` name). An optional second column (tab-separated) gives the RAM (in MB) STAR needs for that pool; empty means the default. Lines starting with `#` and empty lines are skipped.

```text
Sample_P1-cDNA
Sample_P2-cDNA	64000
# Sample_P3-cDNA	(skipped)
```

The fastqs are looked up in `experiment: inputs: fastqs_dir` (per modality) with the aligner's `fold_struct` and `R1_suffix` / `R2_suffix`. See [Wildcard Processing](sub_workflows.md#wildcard-processing) for the YAML and directory forms.

## Wet-lab sheet

`experiment: inputs: wet_lab_info`. Read by the HTO stages: `kite` (the feature barcodes of each pool) and `hashsolo` demultiplexing (each hashtag's sample). Not needed for runs without HTOs.

Tab-separated, one row per sample (hashtag) per pool. The column names are yours; the config says which column is which:

- `modules: kite: columns_to_pick`: pool, hashtag, hashtag barcode sequence
- `modules: hashsolo: columns_to_pick`: pool, hashtag, hashtag barcode sequence, sample ID

```text
unique_sample_ID    hashtag   ab_barcode      SubID
P1                  HTO1      GTCAACTCTTTAGCG  S001
P1                  HTO2      TGATGGCCTATTGGG  S002
```

- Pool names are matched case-insensitively.
- With two hashtags per sample, set `modules: hashsolo: hto_sep` to the separator in the hashtag and barcode columns, e.g. `HTO1_HTO2` and `AATCG..._ATCGG...`.
- `helper_py_scripts/create_wet_lab_info.py` can build this sheet from the wet-lab's files.

## Pool info (vcf_info)

`modules: gt_demux: vcf_info`. Read by genotype demultiplexing (`experiment: demultiplexing: method: vireo` or `both`). vireo always needs it, for the number of donors in each pool. Also read by STARsolo's WASP mode (`modules: STARsolo: wasp: true`), which needs its VCF column.

The columns are found by **position**; their names in the header row don't matter. Each pool appears once, matching the `{pool}` names (case-insensitively for the donor count).

**2 columns:** pool and number of donors. vireo runs without per-pool genotypes; a genotyped run gives vireo the [all-donor VCF](#all-donor-vcf) instead.

```text
pool    n_donors
P1      6
P2      8
```

**4 columns:** pool, number of donors, donors (comma-separated, as in the VCF), and the pool's VCF (bgzipped and indexed).

```text
pool    n_donors   donors                     vcf
P1      3          D1,D2,D3                   /data/vcf/P1.vcf.gz
P2      2          D4,D5                      /data/vcf/P2.vcf.gz
```

- **What the VCF is used for:** `modules: gt_demux: pool_vcf_use` (default `[pileup, vireo]`) sets it. `pileup` makes it the sites of cellSNP's pileup, and `vireo` passes it to vireo as the donor genotypes.
- **STARsolo's WASP mode** needs the 4-column form: every pool aligned needs a VCF.
- **Separator:** `.tsv` is tab-separated, `.csv` comma-separated and `.txt` space-separated.

## All-donor VCF

`modules: gt_demux: donors_per_pool_file`. A single bgzipped, indexed VCF with the genotypes of all donors. vireo reads it when the run is genotyped and the [pool info](#pool-info-vcf_info) has 2 columns. `null` otherwise.

## Reference SNP sites

`modules: gt_demux: genome_1k_ref`. A bgzipped VCF of common SNP sites, e.g. 1000 Genomes with MAF ≥ 5%, on the same chromosome names as the bams. cellSNP piles up these sites when there is no per-pool VCF to use: runs without genotypes, a 2-column pool info file, or `pool_vcf_use` without `pileup`.

## Vireo donor-name converter

`modules: gt_demux: donorName_conv`. Renames vireo's donors (the VCF sample IDs) in the demultiplexed h5ad, e.g. to the study's subject IDs. `file: null` keeps vireo's names.

`.csv` or `.tsv` (or tab-separated `.txt`), with:

| key | |
|---|---|
| `file` | the table |
| `header_lev` | number of header rows (default 1) |
| `pool_col` | column with the pool names (as the run's `{pool}`) |
| `donor_col` | column with vireo's donor names |
| `convert_col` | column with the new names |
| `new_h5ad_colname` | obs column of the h5ad that gets the new names |

```text
pool,vcf_id,SubID
P1,GENO_0001,S001
P1,GENO_0002,S002
```

## Swap correction file

`modules: demultiplex: swap_correction`. A table of the confirmed donors of each pool, e.g. after checking the genotypes (`identify_swaps`). It is read for:

- **the final h5ad of each pool** (`demux_final`): which demultiplexing version to use, and the renamed or dropped donors
- **splitting by those h5ads** (`split_bams: split_by: input: final`)
- **`workflow: start: split_bams`**: its pools are the pools of the run

`.csv` or `.tsv`, one row per pool and donor. The config names the columns:

| key | column |
|---|---|
| `file` | the table |
| `pool_column` | pool |
| `version_column` | the pool's demultiplexing version: `<tool>:<run_tag>` (tool `solo`, `vireo` or `both`, e.g. `vireo:w_gt`; two tools: `solo:hto,vireo:w_gt`), or a name in [`versions`](#demultiplex-versions). One version per pool |
| `donor_from` | the donor as demultiplexed (e.g. the VCF ID) |
| `donor_to` | its confirmed name; empty drops the donor (`Removed`) |
| `status_column` | optional: values other than `True` / `yes` / `1` drop the donor; pools with no donor kept aren't split |

```text
Sample   Donor_name  Solution  final_status  Version_used
P1       donor2      donA      True          w_gt
P1       donor3      donB      True          w_gt
P2       donor4                False         w_gt
```

The stand-alone `helper_py_scripts/donor_bams.py` reads the same kind of table to name per-donor bams. See [Per-donor bams named by the swap correction file](sub_workflows.md#per-donor-bams-named-by-the-swap-correction-file-stand-alone).

## Demultiplex versions

`modules: demultiplex: swap_correction: versions`. What the names in `version_column` mean, when they aren't `<tool>:<run_tag>`:

- **a JSON file** whose `demultiplex_paths` gives each name's `demux` (h5ad) dir and, optionally, its `info` (stats) dir, e.g. `demul_samples_annotate_h5ad_template.json`:

  ```json
  {"demultiplex_paths": {
      "w_gt": {"demux": "/data/demux/vireo_w_gt/", "info": "/data/info/vireo_w_gt/"},
      "wo_gt": {"demux": "/data/demux/vireo_wo_gt/"}}}
  ```

- **a map in the config**, from each name to `<tool>:<run_tag>` (a run of this pipeline) or to `{demux, info, tool}`:

  ```yaml
  versions:
    w_gt: "vireo:w_gt"
    old: {demux: /data/old_demux/, info: /data/old_info/, tool: vireo}
  ```

In each version's dirs, the h5ad and stats file names follow the tool's `fold_struct` and `final_count_matrix_h5ad`. The tool is `split_bams: split_by: demux` unless given.

## External donor labels

`modules: split_bams: split_by: external: file`, with `split_by: input: external`. Donor labels made outside the pipeline, used to split the bams instead of this run's demultiplexing: e.g. a pooled, annotated h5ad of many pools.

A pooled `.h5ad`, or a table (`.csv`, or tab-separated `.tsv` / `.txt`; optionally `.gz`), with one row per cell:

| key | |
|---|---|
| `donor_column` | obs or table column with each cell's donor |
| `pool_column` | obs or table column with each cell's pool of origin |
| `barcode_column` | column with the barcodes; `null`: the obs names (h5ad only) |
| `barcode_transformations` | steps from the file's barcodes to the bam's CB tags |
| `pool_transformations` | steps from the file's pool names to the `{pool}` names |
| `min_barcode_fraction` | least fraction of a bam's sampled reads with a CB tag that are from the labelled cells (default `0.02`); below it the run stops (`check_external_barcodes`) |

```text
barcode                     individualID   pool_of_origin
NPSAD-P12_L1_AAACCTGAGC-1   R1234567       P12
```

Doublets, negatives and unassigned cells are dropped. The run stops when the transformed barcodes don't look like the bam's CB tags. See [Splitting bams by external donor labels](sub_workflows.md#splitting-bams-by-external-donor-labels) for the transformation steps.

## Pools table (bams.tsv)

`modules: split_bams: split_by: external: pools`, with `workflow: start: split_bams`. The bam of each pool, when the bams aren't where the aligner's `bams_dir` and `fold_struct` put them. `null`: the aligner's folders.

Tab-separated, or comma-separated for `.csv`:

| column | |
|---|---|
| `pool` | the `{pool}` name: the labels' pools after `pool_transformations`, or the swap correction file's pools. No `/` |
| `bam` | the pool's coordinate-sorted bam |
| `bai` | optional; default `<bam>.bai` |
| `atac_bam` | multiome (`cellranger_arc`) with `split_bams: modalities: [cDNA, ATAC]`: the pool's ATAC bam (`bam` is then the gex bam) |
| `atac_bai` | optional; default `<atac_bam>.bai` |

```text
pool    bam
samp1   /abs/root/dayA/samp1/samp1_Aligned.sortedByCoord.out.bam
samp3   /abs/root/dayB/samp3/samp3_Aligned.sortedByCoord.out.bam
```

- **Pools of the labels that aren't in the table** aren't split, with a warning. If none of them is in the table, the run stops.
- **A missing bam or bai** stops the run, with every missing file listed.
- **The bams are only read.** Derived bams go to `split_bams: temp_dir`.
- **Multiome** (`cellranger_arc`): `bam` is the gex bam; with ATAC in `split_bams: modalities`, every pool needs an `atac_bam` too (a missing one stops the run).

See [Bams outside the aligner's folders](sub_workflows.md#bams-outside-the-aligners-folders-pools-table) for a one-liner that builds the table.

## Split-bam donor-name converter

`modules: split_bams: donor_name_converter`. Renames the donors before the bams are split, so that each per-donor bam is named by the new name. `file: null` keeps the demultiplexed names.

A **comma-separated** table (read as CSV, whatever its extension). `from_column` holds the donor names as demultiplexed or labelled, and `to_column` the new names. With external labels, donors missing from the table keep their names.

```text
vcf_id,individualID
GENO_0001,R1234567
GENO_0002,R7654321
```

## Reference genotypes for MBV

`modules: identify_swaps: ref_vcf`. A bgzipped, indexed VCF with the genotypes of the expected donors. QTLtools mbv checks each per-donor bam against all of them, which shows sample swaps. The bams are the run's split bams (`mbv_inp: vireo_outs`); otherwise, `mbv_inp` is a path template of other bams.

## WASP samples table (wasp_samples.tsv)

`modules: wasp: samples`. Required for `workflow: until: wasp` and `start: wasp`. It lists the per-donor bams to WASP-filter, and each donor's genotype ID. Tab-separated, or comma-separated for `.csv`, with one row per donor bam (`helper_py_scripts/make_wasp_samples.py` writes it from per-donor bams, e.g. split bams: see [Helper scripts](helper_scripts.md)):

| column | required | |
|---|---|---|
| `pool` | always | pool (no `/`) |
| `donor` | always | donor (no `/`) as the final demultiplexing named it: after `split_bams`, the split bam's donor (e.g. its VCF ID) |
| `final_donor` | no | the donor after swap correction (no `/`); names the output bam in its pool's folder (`<wasp_dir><modality>/<pool>/<final_donor>.after_wasp.bam`), once per pool and modality. Swap corrections are only applied to the h5ad files so far, so the split bams keep the demultiplexed name. Empty, or no column: from the [sample names](#wasp-sample-names), else the donor. With no `final_donor` column in either table, a warning says no swap corrections are assumed |
| `modality` | no | `cDNA` or `ATAC` (any case). Empty, or no column: with `start: wasp`, `cDNA`; after `split_bams`, every split modality of the donor (`split_bams: modalities`), unless the donor has a row of its own for that modality |
| `vcfID` | no | the final donor's genotype sample ID(s) in its genotype set (as in the VCF; a `sample` column is read too), separated by `;` (or `,` in a tab-separated table), with no whitespace in an ID. They're written one per line to the file given as `--samples` to `find_intersecting_snps.py` (WASP takes a file, whose lines it reads up to the first whitespace, or a comma-separated string). Empty, or no column: from the [sample names](#wasp-sample-names), else the final donor's name |
| `bam` | `start: wasp` | the per-donor bam (coordinate-sorted) |
| `bai` | no | default `<bam>.bai` |
| `genotype_set` | with more than one set | the [genotype set](#wasp-genotype-sets) holding `vcfID` |

```text
pool    donor   final_donor  vcfID       bam                                  genotype_set
P1      D1      D1           GENO_0001   /data/split/Sample_P1-cDNA/D1.bam    wgs
P1      D2      D3           GENO_0003   /data/split/Sample_P1-cDNA/D2.bam    array
```

Multiome samples (`start: wasp`): one row per bam, e.g. the gex and ATAC bams of a cellranger-arc run:

```text
pool    donor   final_donor  vcfID       modality  bam
S001    S001    S001         NDA_0001    cDNA      /data/arc/S001/outs/gex_possorted_bam.bam
S001    S001    S001         NDA_0001    ATAC      /data/arc/S001/outs/atac_possorted_bam.bam
```

**Checks:**

- Missing columns, an empty pool or donor, `/` in names, a modality other than `cDNA` or `ATAC`, and duplicated (pool, donor, modality) rows stop the run.
- So do two donors of a pool and modality with the same `final_donor` (an empty one counts as the donor).
- Both a `vcfID` and a `sample` column stop the run (keep `vcfID`).
- So does a `genotype_set` that isn't in `modules: wasp: genotypes`.
- Every `sample` must be in its set's genotypes. A rule checks this once, before any WASP job.
- After `split_bams`, split donors not in the table aren't filtered, with a warning. If none is in the table, the run stops.
- With `start: wasp`, a missing bam or bai stops the run, with all of them listed.

## WASP sample names

`modules: wasp: sample_names`. Optional. The final donor and genotype sample ID(s) of each demultiplexed donor (or of every donor of a pool). The bams are split by the final demultiplexing, but swap corrections are only applied to the h5ad files so far. So a donor whose labels were swapped needs its corrected name, which names its WASP outputs, and the VCF ID of that corrected donor. It also covers pools or donors not named as in the genotypes, e.g. one donor per pool, named after the pool. It fills the [samples table](#wasp-samples-table-wasp_samplestsv)'s rows with no `final_donor` or `vcfID`. Tab-separated, or comma-separated for `.csv`:

| column | required | |
|---|---|---|
| `pool` | always | pool |
| `donor` | no | the donor, as the final demultiplexing named it; empty, or no column: every donor of the pool |
| `final_donor` | no | the donor after swap correction (no `/`); empty: the donor |
| `vcfID` | always | genotype sample ID(s) of the final donor, as in the VCF (`;`- or `,`-separated, no whitespace). A `sample` column, as in tables from before `vcfID`, is read too |

```text
pool       donor   final_donor  vcfID
Pool_A     D1      D1           SUBJ_0101
Pool_A     D2      D3           SUBJ_0103
Pool_A     D3      D2           SUBJ_0102
Pool_B             Pool_B       SUBJ_0201
```

Here `Pool_A`'s `D2` and `D3` were swapped: `D2`'s split bam is WASP-filtered with `SUBJ_0103`'s genotypes and written as `<wasp_dir>cDNA/Pool_A/D3.after_wasp.bam`. As with the split bams, each pool has its own folder, so a final donor in several pools (e.g. a donor run in two pools) gets one bam per pool, none overwriting another; within a pool and modality a final donor can only be given once.

**No swap corrections:** without a `final_donor` column, here or in the samples table, a warning says that no swap corrections are assumed, and `final_donor` is copied from `donor`. Outputs are then named as before, so a run looks the same whether or not swaps were corrected:

```text
Warning: WASP samples table wasp_samples.tsv and sample names table names.tsv have no final_donor column: no swap corrections assumed (final_donor = donor)
```

Rows whose `final_donor` is empty while the column exists are listed in a warning too.

A row's `final_donor` and `vcfID`, in order: the samples table's; else the (pool, donor) row here; else the pool's row (empty donor); else the donor (for `vcfID`: the final donor's name). An empty `pool` or `vcfID`, `/` in `final_donor` and duplicated (pool, donor) rows stop the run. The IDs are checked against the genotypes like the table's.

## WASP genotype sets

`modules: wasp: genotypes`. Named sets, e.g. one for WGS and one for SNP arrays. Each set is one of:

- **`h5_dir`:** a folder with WASP's `snp_index.h5`, `snp_tab.h5` and `haps.h5`, made with WASP's `snp2h5`.
- **`vcf` and `chrom_sizes`:** an indexed VCF/BCF, or a list of them (e.g. one per chromosome). `chrom_sizes` is a two-column file of chromosome and length, as for `snp2h5 --chrom`. Rule `wasp_snp2h5` converts the set once to `<wasp_dir>genotypes/<set>/`, keeping the set's samples in the table and the sites variable among them.

```yaml
modules:
  wasp:
    genotypes:
      wgs: {h5_dir: /data/wasp_h5/wgs/}
      array: {vcf: /data/array.bcf, chrom_sizes: /data/hg38.chrom.sizes}
```

**Chromosome names:** the genotypes and the bams must name chromosomes the same way (`chr1` or `1`). Each job checks this before it starts.

## WASP remap references

The references the reads over SNPs are remapped to. They must be those the bams were aligned to.

- **cDNA:** `modules: wasp: rna: remap`: the STAR index (`genome`) and, optionally, the GTF (`gtf`). After `split_bams`, `null` takes the aligner's reference: STARsolo's `genome_pick`, or the `star` folder of the cellranger / cellranger-arc reference; the latter may need a STAR index rebuilt for a newer STAR (see [config caveats](config_caveats.md)).
- **ATAC:** `modules: wasp: atac: remap: bwa_index`: the bwa-indexed fasta. After `split_bams`, `null` takes the cellranger-arc reference's `fasta/genome.fa` (`cellranger-arc mkref` indexes it for bwa).
- **With `start: wasp`:** each modality of the samples table needs its own (`genome` for cDNA rows, `bwa_index` for ATAC rows).

See [WASP filtering of per-donor bams](sub_workflows.md#wasp-filtering-of-per-donor-bams) for the jobs, outputs and usage.
