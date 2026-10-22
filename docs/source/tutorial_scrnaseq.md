# scRNA-seq (one wildcard)

This tutorial aligns samples of one donor each (no demultiplexing), collects their QC metrics, and writes an h5ad per sample. Each sample is one value of the `pool` wildcard; see {doc}`tutorial_pooled_snrnaseq` for the inputs (fastqs and `select_fastqs`).

## Aligner

The prefix of `last_step` picks the aligner:

- `STARsolo_*`: STARsolo
- `cellranger_*`: cellranger count (`cellranger_count` in the config)
- `multiome_*`: cellranger-arc count, for multiome runs (cDNA and ATAC)

## Alignment and QC metrics

```yaml
last_step: STARsolo
picard_metrics: all    # all, gc, rnaseq or null
STARsolo_pipeline:
  count_matrix: Gene   # Gene (exons, e.g. cells) or GeneFull (exons and introns, e.g. nuclei)
```

With `STARsolo_pipeline: features` listing both (`Gene GeneFull …`), STARsolo writes both count matrices; `count_matrix` picks the one used after the alignment. Picard's CollectGcBiasMetrics and CollectRnaSeqMetrics run on each bam (STARsolo runs only).

## h5ad of each sample

```yaml
last_step: STARsolo_demux
demux_type: null
create_final_h5ad: True
max_mito_percentage: 5
min_genes_per_cell: 1000
min_cells_per_gene: 10
```

Without demultiplexing, `demux_run_none` writes the h5ad of each sample with the QC filters (cells and genes that pass them are marked in `QC_pass`) and `demux_final` writes it to `demultiplex: demux_count_matrix_dir`, with the QC stats in `demultiplex: demultiplex_info_dir`.

With cellranger count instead:

```yaml
last_step: cellranger_demux
demux_type: null
cellranger_count:
  bams_dir: projdir/cellranger_count_outputs/
  ref: *cellR_gex2024
  include_introns: false   # true for nuclei
```

## Allele-specific mapping bias (WASP)

For analyses of allele-specific expression, STAR's WASP mode flags reads with mapping bias (SAM attribute `vW`). It needs the sample's vcf (the `vcf` column of `gt_demux_pipeline: vcf_info`):

```yaml
STARsolo_pipeline:
  wasp: true
gt_demux_pipeline:
  vcf_info: projdir/sample_info.tsv
```

## Collecting the stats

`helper_sh_scripts/run_update_logs.sh` compiles STARsolo's and Picard's stats of all samples, including STARsolo's run time (`RUN_TIME`), into one table.
