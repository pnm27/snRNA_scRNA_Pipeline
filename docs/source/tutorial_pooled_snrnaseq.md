# Pooled snRNA-seq (one wildcard)

This tutorial demultiplexes pools of nuclei from several donors with their genotypes (cellSNP and vireoSNP), checks the donors for swaps, and writes one bam per donor. Each pool is one value of the `pool` wildcard.

It runs the pipeline in rounds, each with its own `run_tag`, so the outputs of earlier rounds are kept:

1. demultiplex without genotypes (`wo_gt`)
2. check the donors of each pool against their genotypes (identify swaps)
3. demultiplex with genotypes (`w_gt`)
4. pick the final demultiplexing of each pool, correct swapped donors, and split the bams by donor

## Inputs

Fastqs of each pool are in their own folder; the folder and file names carry the pool name (here `P1`), followed by `-cDNA` (`-HTO` for hashtag fastqs):

```text
projdir/fastq_files/
├── Sample_P1-cDNA/
│   ├── P1-cDNA_S1_L001_R1_001.fastq.gz
│   └── P1-cDNA_S1_L001_R2_001.fastq.gz
└── Sample_P2-cDNA/
    └── ...
```

`select_fastqs` (here `projdir/fastq_files.txt`) lists the pools, one per line. An optional second column (tab-separated) has the memory (MB) for STARsolo of that pool:

```text
P1-cDNA
P2-cDNA	60000
```

`gt_demux_pipeline: vcf_info` describes the donors of each pool. It has headers; the column names are free, but their order is fixed:

- without genotypes: `pool, n_donors`
- with genotypes: `pool, n_donors, donors, vcf` (`vcf`: the pool's donors' genotypes)

```text
SetNumber	Gt_count	donors	pool_vcf
P1	4	D1,D2,D3,D4	projdir/vcfs/P1.vcf.gz
```

## Configuration

Copy `new_config.example.yaml` to `new_config.yaml` and set the project's paths (`select_fastqs`, `cDNA_fastqs_dir`, the `*_dir` keys), the genome (`STARsolo_pipeline: genome_pick`) and the folder structures (`fold_struct*`, see {doc}`sub_workflows`). The config is checked against `schema/config.schema.yaml` at the start of every run.

Check what a run will do with a dry run before running it:

```bash
snakemake -n
```

## Round 1: without genotypes

```yaml
last_step: STARsolo_gt_demux
demux_type: vireo
demux_run_type: wo_gt
run_tag: wo_gt
create_final_h5ad: True
```

STARsolo aligns each pool, cellSNP piles up the common SNPs of `genome_1k_ref` in its cells, and vireo assigns the cells to `n_donors` donors (`donor0`, `donor1`, …). As `run_tag` is `wo_gt`, the cellSNP, vireoSNP, count matrix and demultiplex info dirs each get a `wo_gt/` sub-folder.

The outputs of each pool:

- `final_count_matrix/vireoSNP/wo_gt/…h5ad`: the h5ad of this run (QC stats in `uns['qc_stats']`), and `demultiplex/info/vireoSNP/wo_gt/…_vireo_stats.tsv`
- `final_count_matrix/final/<pool>_final.h5ad`: the final h5ad, with the donor of each cell in `Donor_final`
- `demultiplex/info/final/`: the basic stats and demultiplex info of each pool, and `final_demux_versions.tsv` (one row per pool)

## Round 2: identify swaps

```yaml
last_step: STARsolo_gt_demux_identify_swaps
gt_check: yes
split_bams_pipeline:
  subset_chr: 1
  split_by:
    demux: vireo
    input: raw
```

Keep `run_tag: wo_gt`, so the bams are split by round 1's vireo output. The bam of each pool is subset to chromosome 1 and split into one bam per donor (`split_bams_dir/gt_check/`), and QTLtools-mbv compares each donor's reads with the genotypes in `identify_swaps: ref_vcf`. Use its outputs (`identify_swaps: mbv_out_dir`) to match `donor0`, `donor1`, … to the genotyped donors, and to find swapped or missing donors.

## Round 3: with genotypes

```yaml
last_step: STARsolo_gt_demux
demux_type: vireo
demux_run_type: w_gt
run_tag: w_gt
gt_demux_pipeline:
  vcf_info: projdir/pool_info.tsv    # with the vcf column
  pool_vcf_use: [pileup, vireo]
  reuse_cellsnp_dir: null
```

vireo now assigns the cells to the genotyped donors. `pool_vcf_use` sets what the pool's vcf is used for: the cellSNP pileup and/or vireo's donor genotypes. If the SNPs didn't change since an earlier run, `reuse_cellsnp_dir` (that run's `cellsnp_dir`, e.g. `projdir/demultiplex/cellSNP/wo_gt/`) runs vireo on its cellSNP outputs.

## Round 4: final h5ad and per-donor bams

The swap correction file lists, for each pool, the demultiplexing run to use and the donors to rename (one row per pool and donor). An empty `donor_to` drops that donor (its cells' `Donor_final` is `Removed`):

```text
pool	version	donor_from	donor_to
P1	vireo:w_gt	D2	D3
P1	vireo:w_gt	D3	D2
P2	vireo:wo_gt	donor1
```

`version` is `<tool>:<run_tag>`, with tool `vireo`, `solo` or `both` (both tools of that run); runs of different tools are separated by commas, e.g. `solo:hto,vireo:w_gt`.

```yaml
last_step: STARsolo_split_bams_gt_demux
gt_check: no
demultiplex:
  swap_correction:
    file: projdir/swap_correction.tsv
    pool_column: pool
    version_column: version
    donor_from: donor_from
    donor_to: donor_to
split_bams_pipeline:
  split_by:
    demux: vireo
    input: final
```

`demux_final` rebuilds the final h5ad of each pool from the run in the swap correction file, and the bam of each pool, filtered by cell barcode, is split by `Donor_final` into `split_bams_dir/final/`.

## Adding hashtag (HTO) demultiplexing

With hashtag fastqs (`HTO_fastqs_dir`) and the wet lab file (`wet_lab_info`), hashsolo can run with vireo (`demux_type: both`), or be added to an earlier vireo run:

```yaml
demux_type: add_solo
run_tag: hto
demultiplex:
  previous_run_tag: w_gt
```

`demux_final` then combines both: hashsolo's donors (`SubID_cs`) and vireo's (`SubID_vs`) in the same h5ad (`Donor_final` follows vireo). Pools in the swap correction file use the runs listed there, e.g. `solo:hto,vireo:w_gt`.

## Collecting the stats

`helper_sh_scripts/run_update_logs.sh` compiles STARsolo's, Picard's and the demultiplexing stats of each pool into one table, including STARsolo's run time (`RUN_TIME`). The demultiplex info of each final h5ad (`demultiplex/info/final/<pool>_final_info.tsv`) has the basic and tool stats together, in the format it reads.
