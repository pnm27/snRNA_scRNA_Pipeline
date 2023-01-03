# snRNA_scRNA_Pipeline Introduction
Snakemake Pipeline with "selectable" modules for snRNA seq pre-processing. It also supports various software/pipeline for scRNA seq pre-processing.

The highlights of the pipeline are:
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
	<li><a href="https://scanpy.readthedocs.io/en/stable/" name="sc">Scanpy Manual</a></li>
	<li><a href="https://snakemake.readthedocs.io/en/stable/" name="snk">Snakemake Manual</a></li>
	<li><a href="https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md" name="sts">STARsolo Manual</a> </li>
	<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics" name="gcb">GC bias metrics (PICARD)</a></li>
	<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics" name="rna">RNA seq metrics (PICARD)</a></li>
	<li><a href="https://github.com/pachterlab/kite"name="kite">KITE <i>(kallisto indexing and tag extraction)</i></a></li>
	<li><a href="https://cellsnp-lite.readthedocs.io/en/latest/manual.html" name="csp">cellSNP Manual</a></li>
	<li><a href="https://vireosnp.readthedocs.io/en/latest/manual.html" name="cir">vireoSNP Manual</a></li>
	<li><a href="https://github.com/calico/solo#how-to-demultiplex-cell-hashing-data-using-hashsolo-cli" name="hsolo">hashsolo Info</a></li>
	<li><a href="https://qtltools.github.io/qtltools/pages/QTLtools-mbv.1.html" name="mbv">QTLtools-mbv</a></li>
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

```sh
sh run_snakemake.sh
```

## Configuration for pipeline
This pipeline depends on a yaml config file (new_config.yaml), which has all relevant options for each rule present in this pipeline. This not only streamlines the process of maintaining or modifying program specific parameters but also makes repetitive usage (can be within the project or using it over mutliple projects) of this pipeline over time very efficient.