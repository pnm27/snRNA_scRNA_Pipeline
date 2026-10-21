# Overview of the Pipeline
To set up the pipeline for projects, the following setups need to be performed:

1. Setting up wildcards for the project (see {ref}`how <sub_workflows:Wildcard Processing>`).
2. According to the experiment select a module (module options).
3. As required and wherever possible, setup the folder structures for different programs.
4. Check and modify parameters of the programs to be run (in new_config.yaml file).
   1. (Currently) make sure all required python and R packages are present.
5. Change global variables as required (in Snakefile - check).

# Wildcard Processing
For the purpose of creating wildcards a list of samples to be processed is provided to the pipeline. There are 3 ways to achieve this:

- list of samples/pools (as a folder structure)
- yaml file containing the list of samples/pools
- Directory containing input files

## List of samples



``` {currentmodule} helper_py_scripts/

```

This pipeline has many combinations of the aforementioned programs as a built-in set that can be executed using specific keywords.

```{eval-rst}
.. _selectable-modules:
```

# Selectable Modules
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
* - all_multi_vcf
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_kb_solo
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
* - demultiplex.smk
  - It contains the rules that create an h5ad per demultiplexing tool (`demux_run_<tool>`) and the final h5ad of each pool (`demux_final`). See [Demultiplexing outputs](#demultiplexing-outputs).
* - identify_swaps.snkmk 
  - It contains rules for identifying swaps using QTLtools-mbv.
```

[^consol]: Consolidating into one

[^not-yet]: Not yet implemented

# Demultiplexing outputs

Each demultiplexing run writes to its own sub-folder (`run_tag` in the config, e.g. `wo_gt` or `w_gt_2`) of the calico_solo, cellSNP, vireoSNP, final count matrix and demultiplex info dirs. With `run_tag: null` the dirs are used as they are.

```{list-table} Demultiplexing rules
:header-rows: 1
:name: demux-rules-table

* - Rule
  - Outputs
* - `demux_run_solo`, `demux_run_vireo` (`demux_run_none` without demultiplexing)
  - An h5ad with the classification of one tool (QC stats in `uns['qc_stats']`) and the stats of that tool (`_hashsolo_stats.tsv` / `_vireo_stats.tsv`), in the tool's dirs.
* - `demux_final`
  - The final h5ad of each pool (`demultiplex: demux_count_matrix_dir`, `fold_struct_final`), with the donor of each cell in `Donor_final`; its basic stats (`_basic_stats.tsv`: QC and the demultiplex versions used) and its demultiplex info (basic and tool stats together, read by `update_logs.py`).
* - `collect_demux_versions`
  - `final_demux_versions.tsv` in `demultiplex: demultiplex_info_dir`: one row per pool with the demultiplex versions used, whether the swap correction was applied, and the number of donors and of cells passing QC.
```

`demux_final` combines, for each pool:

1. the runs in its row of the swap correction file (`demultiplex: swap_correction`), given as `<tool>:<run_tag>` (tool: `solo`, `vireo` or `both`), whose `donor_from` → `donor_to` columns rename the donors (an empty `donor_to` marks the donor's cells `Removed`), else
2. the tools of this run and, for `add_solo` / `add_vireo`, the other tool from the run `demultiplex: previous_run_tag`.

To split the bams by the final h5ad (e.g. after the swap correction), set `split_bams_pipeline: split_by: input: final`.

Split bams go to a sub-folder of `inp_split_bams_dir`, `split_bams_dir`, `per_donor_split_log_dir` and `new_temp_dir`: `gt_check/` for genotype checks (`gt_check: yes`; the whole bam, or the one subset to `subset_chr`) and `final/` otherwise (the bam filtered by cell barcode).

```{eval-rst}
.. autosummary::
  :nosignatures:

   demultiplex_helper_funcs._ret_htos_calico_solo  
```
