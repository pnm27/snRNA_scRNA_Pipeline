# scRNA-seq (one wildcard)

This tutorial aligns samples of one donor each (no demultiplexing), collects their QC metrics, and writes an h5ad per sample. Each sample is one value of the `pool` wildcard; see {doc}`tutorial_pooled_snrnaseq` for the inputs (fastqs and `select_fastqs`).

## Aligner

`workflow: aligner` picks the aligner, and `workflow: until` the last stage to run (see {doc}`sub_workflows`):

- `STARsolo`: STARsolo
- `cellranger`: cellranger count (`modules: cellranger` in the config)
- `cellranger_arc`: cellranger-arc count, for multiome runs (cDNA and ATAC)

## Alignment and QC metrics

```yaml
workflow:
  aligner: STARsolo
  until: align
  picard_metrics: all    # all, gc, rnaseq or null
modules:
  STARsolo:
    count_matrix: Gene   # Gene (exons, e.g. cells) or GeneFull (exons and introns, e.g. nuclei)
```

With `modules: STARsolo: features` listing both (`Gene GeneFull …`), STARsolo writes both count matrices; `count_matrix` picks the one used after the alignment. Picard's CollectGcBiasMetrics and CollectRnaSeqMetrics run on each bam (STARsolo runs only).

## h5ad of each sample

```yaml
workflow:
  until: demux
experiment:
  demultiplexing:
    method: null
    create_final_h5ad: true
  qc:
    max_mito_percentage: 5
    min_genes_per_cell: 1000
    min_cells_per_gene: 10
```

Without demultiplexing, `demux_run_none` writes the h5ad of each sample with the QC filters (cells and genes that pass them are marked in `QC_pass`) and `demux_final` writes it to `modules: demultiplex: demux_count_matrix_dir`, with the QC stats in `modules: demultiplex: demultiplex_info_dir`.

With cellranger count instead:

```yaml
workflow:
  aligner: cellranger
  until: demux
modules:
  cellranger:
    bams_dir: projdir/cellranger_count_outputs/
    ref: *cellR_gex2024
    include_introns: false   # true for nuclei
```

## Allele-specific mapping bias (WASP)

For analyses of allele-specific expression, STAR's WASP mode flags reads with mapping bias (SAM attribute `vW`). It needs the sample's vcf (the `vcf` column of `modules: gt_demux: vcf_info`):

```yaml
modules:
  STARsolo:
    wasp: true
  gt_demux:
    vcf_info: projdir/sample_info.tsv
```

## Collecting the stats

`helper_sh_scripts/run_update_logs.sh` compiles STARsolo's and Picard's stats of all samples, including STARsolo's run time (`RUN_TIME`), into one table.
