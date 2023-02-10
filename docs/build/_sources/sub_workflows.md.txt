# Modules


``` {currentmodule} helper_py_scripts/

```

This pipeline has many combinations of the aforementioned programs as a built-in set that can be executed using specific keywords. 

```{eval-rst}
.. _selectable-modules:
```

## Selectable Modules
The following combinations of programs can be run:

<ul>
<li> all</li>
<li> starsolo</li>
<li> starsolo_rnaseqmet</li>
<li> starsolo_gcbiasmet</li>
<li> starsolo_kb_solo</li>
<li> starsolo_picard</li>
<li> starsolo_gt_demux</li>
<li> starsolo_split_bams</li>
<li> starsolo_split_bams_gt_demux</li>
<li> starsolo_split_bams_gt_demux_multi_vcf</li>
<li> starsolo_gt_demux_multi_vcf</li>
<li> starsolo_cellsnp</li>
<li>starsolo_rnaseqmet_kb_solo</li>
<li>starsolo_gcbiasmet_kb_solo</li>
<li>starsolo_gt_demux_identify_swaps</li>
<li>starsolo_resolve_swaps_gt_demux <sup><small>#</small></sup></li>
</ul>

where **starsolo** represents STARsolo; **rnaseqmet** and **gcbiasmet** refer to PICARD's CollectRnaSeqMetrics and CollectGcBiasMetrics, respectively while **picard** represents inclusion of both the previously-mentioned programs; **kb_solo** refers to using kallisto, bustools and calico_solo for demultiplexing; **gt_demux** refers to using cellSNP and vireoSNP for genotype based demultiplexing; **split_bams** refers to splitting pooled/multiplexed bams using hashsolo's outputs while **split_bams_gt_demux** refers to splitting pooled/multiplexed bams using vireo's output; **identify_swaps** refers to using qtltools_mbv. The option **multi_vcf** is to provide muiltiple runs (i.e. multiple sets of vcf inputs) for the same sample.

## Module description

```{list-table} Modules_info
:header-rows: 1
:name: module-info-table

* - Module Name
  - Module Info
  - Sub Worflows Involved
* - all
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_rnaseqmet
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gcbiasmet
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_kb_solo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_picard
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gt_demux
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_split_bams
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_split_bams_gt_demux
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_split_bams_gt_demux_multi_vcf
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gt_demux_multi_vcf
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_cellsnp
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_rnaseqmet_kb_solo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gcbiasmet_kb_solo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gt_demux_identify_swaps
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_resolve_swaps_gt_demux
  - module_info (more desc in its own file)
  - sub_wkfl
```

# Sub-Snakemake workflows
This pipeline divides each module into its self-contained individual workflows. These are:

```{list-table} Sub_Snakemake_workflows_table
:header-rows: 1
:name: sub_smk-info-table

* - Name of Workflow
  - Description
* - resources.snkmk 
  - It contains memory (in MB per thread) and time requirements (in minutes) for each rule.
* - calico_solo_demux.snkmk 
  - It contains hashsolo rule.
* - split_bams.snkmk [^consol]
  - It contains rules needed to split pooled bams into individual bams dependent on output
  produced by either hashsolo or vireoSNP using custom scripts.
* - input_processing.snkmk 
  - It contains rules that collects values for all the wildcards.
* - STARsolo.snkmk 
  - It contains rules for STARsolo.
* - produce_targets.snkmk 
  - It contains the **rule all** and the needed functions.
* - snv_aware_align.snkmk [^not-yet]
  - **This might be removed soon**
* - kite.snkmk  
  - It contains rules for the [kite](https://github.com/pachterlab/kite#kite-kallisto-indexing-and-tag-extraction) workflow.
* - picard_metrics.snkmk 
  - It contains rules for all PICARD metrics (GCBiasMetrics and RNAseqMetrics).
* - pheno_demux3.snkmk 
  - It contains rules for the cellSNP-vireoSNP pipeline.
* - split_bams_gt.snkmk [^consol]
  - It contains rules needed to split pooled bams into individual bams dependent on output produced by vireo.
* - demultiplex_no_argp.snkmk
  - It contains rules for demultiplexing using hashsolo and/or vireoSNP output and create a count matrix file.
* - identify_swaps.snkmk 
  - It contains rules for identifying swaps using QTLtools-mbv.
```

[^consol]: Consolidating into one

[^not-yet]: Not yet implemented

```{eval-rst}
.. autosummary::  
  :toctree: generated/
  :nosignatures:

   demul_samples_no_argp.ret_htos_calico_solo  
```

## Overview of the pipeline

Directed Acyclic Graph of the whole pipeline:
![DAG](../../images/Whole_pipeline.png)

