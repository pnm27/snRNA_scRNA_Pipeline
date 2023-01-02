# scRNA_snRNA_Pipeline
Snakemake Pipeline with "selectable" modules for scRNA/snRNA seq pre-processing.
The highlights of the pipeline is:
<ul>
	<li> Streamlined processes to modify parameters for each program through a single yaml file </li>
	<li> Easily modifiable to accomodate more rules </li>
	<li> Can be used for both individual samples as well as multiplexed pools </li>
	<li> Preserve folder structures (mirroring fastqs' folder structures) </li>
	<li> Organize outputs from each module </li>
	<li> Select multiple pre-set modules that simplifies usage across multiple projects </li>
</ul>

## Requirements
This pipeline depends on the following packages/programs:
<ul>
	<li><a href="https://scanpy.readthedocs.io/en/stable/">Scanpy Manual</a></li>
	<li><a href="https://snakemake.readthedocs.io/en/stable/">Snakemake Manual</a></li>
	<li><a href="https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md">STARsolo Manual</a> </li>
	<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics">GC bias metrics (PICARD)</a></li>
	<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics">RNA seq metrics (PICARD)</a></li>
	<li><a href="https://github.com/pachterlab/kite">KITE <i>(kallisto indexing and tag extraction)</i></a></li>
	<li><a href="https://cellsnp-lite.readthedocs.io/en/latest/manual.html">cellSNP Manual</a></li>
	<li><a href="https://vireosnp.readthedocs.io/en/latest/manual.html">vireoSNP Manual</a></li>
	<li><a href="https://github.com/calico/solo#how-to-demultiplex-cell-hashing-data-using-hashsolo-cli">hashsolo Info</a></li>
</ul>

### Packages installed through conda
All the packages installed through anaconda3/2018.12 for Python 3.9.5 are described (add_link and file)

### Packages installed through pip
All the packages installed for R version 4.1.0 are described (add_link and file)

## Overview of the pipeline
Directed Acyclic Graph of the whole pipeline:
![DAG](images/Whole_pipeline.png)

## Setting up profiles
The info for setting up profiles for different workload managers is mentioned [here](https://github.com/Snakemake-Profiles)

## Executing Pipeline
This pipeline can be executed by executing (in case of any workflow manager, submitting) the script called `run_snakemake.sh`

```Shell
sh run_snakemake.sh
```

## Configuration for pipeline
This pipeline depends on a yaml config file (new_config.yaml), which has all relevant options for each rule present in this pipeline. This not only streamlines the process of maintaining or modifying program specific parameters but also makes repetitive usage (can be within the project or using it over mutliple projects) of this pipeline over time very efficient.

## Selectable Modules
The following combinations of programs can be run, where starsolo represents STARsolo, rnaseqmet and gcbiasmet refer to PICARD's CollectRnaSeqMetrics and CollectGcBiasMetrics, respectively while picard represents inclusion of both programs, kb_solo refers to using kallisto, bustools and calico_solo for demultiplexing, gt_demux refers to using cellSNP and vireoSNP for genotype based demultiplexing, split_bams_kb_solo refers to splitting pooled/multiplexed bams

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
</ul>

NOTE: These module names are case-insensitive!

### Module Specifications
<dl>
	<dt>All</dt>
	<dd>This module includes alignment through STARsolo and PICARD's both programs (GcBiasMetrics and RnaSeqMetrics) along with calico solo demultiplexing through kallisto bustools output</dd>
	<dt>STARsolo</dt>
	<dd>Only alignemnt through STARsolo.</dd>
	<dt>STARsolo_rnaseqmet</dt>
	<dd>This module will execute STARsolo and PICARD's RNAseq metrics.</dd>
	<dt>STARsolo_gcbiasmet</dt>
	<dd>This module will execute STARsolo and PICARD's GC bias metrics.</dd>
	<dt>STARsolo_kb_solo</dt>
	<dd>This module will execute STARsolo, kallisto_bustools, and demultiplexing by calico_solo/hashsolo.</dd>
	<dt>STARsolo_PICARD</dt>
	<dd>This module will execute STARsolo, PICARD's RNAseq metrics, and PICARD's GC bias metrics.</dd>
	<dt>STARsolo_gt_demux <i>(Not yet implemented)</i></dt>
	<dd>This module will execute cellSNP and vireoSNP</dd>
</dl>


## Sub-Snakemake workflows
This pipeline divides each module into its self-contained individual workflows. These are:

<dl>
	<dt> resources.snkmk </dt>
	<dd> </dd>
	<dt> calico_solo_demux.snkmk </dt>
	<dd> </dd>
	<dt> split_bams.snkmk<sup><a href="#ft1" id="ref1" >*</a></sup></dt>
	<dd> </dd>
	<dt> input_processing.snkmk </dt>
	<dd> </dd>
	<dt> STARsolo.snkmk </dt>
	<dd> </dd>
	<dt> produce_targets.snkmk </dt>
	<dd> </dd>
	<dt> snv_aware_align.snkmk (Not yet implemented)</dt>
	<dd> </dd>
	<dt> kite.snkmk </dt>
	<dd> </dd>
	<dt> picard_metrics.snkmk </dt>
	<dd> </dd>
	<dt> pheno_demux3.snkmk </dt>
	<dd> </dd>
	<dt> split_bams_gt.snkmk<sup><a href="#ft1" id="ref1" >*</a></sup></dt>
	<dd> </dd>
	<dt> demultiplex_no_argp.snkmk </dt>
	<dd> </dd>

<sup id="#ft1">*</sup> Consolidating into one
</dl>

