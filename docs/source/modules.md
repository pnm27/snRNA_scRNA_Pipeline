# Modules
This pipeline has many combinations of the aforementioned programs as a built-in set that can be executed using specific keywords. 

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


# Sub-Snakemake workflows
This pipeline divides each module into its self-contained individual workflows. These are:

<dl>
	<dt> resources.snkmk </dt>
	<dd> This snakemake sub-workflow contains memory (in MB per thread) and time requirements 
	(in minutes) for each rule. </dd>
	<dt> calico_solo_demux.snkmk </dt>
	<dd> This snakemake sub-workflow contains hashsolo rule. </dd>
	<dt> split_bams.snkmk<sup><a href="#ft1" id="ref1" >*</a></sup></dt>
	<dd> This snakemake sub-workflow contains rules needed to split pooled bams into individual
	bams dependent on output produced by hashsolo. </dd>
	<dt> input_processing.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules that collects values for all the wildcards. </dd>
	<dt> STARsolo.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules for STARsolo. </dd>
	<dt> produce_targets.snkmk </dt>
	<dd> This snakemake sub-workflow contains the **rule all** and the needed functions. </dd>
	<dt> snv_aware_align.snkmk <sup><a href="#ft2" id="ref2" >#</a></sup></dt>
	<dd> This snakemake sub-workflow contains rules for </dd>
	<dt> kite.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules for the <a href="https://github.com/pachterlab/kite#kite-kallisto-indexing-and-tag-extraction">kite</a> workflow 
	</dd>
	<dt> picard_metrics.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules for all PICARD metrics mentioned above (<a href="#gcb">GCBiasMetrics</a>
	 and <a href="#rna">RNAseqMetrics</a>)</dd>
	<dt> pheno_demux3.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules for the cellSNP-vireoSNP pipeline. </dd>
	<dt> split_bams_gt.snkmk<sup><a href="#ft1" id="ref1" >*</a></sup></dt>
	<dd> This snakemake sub-workflow contains rules needed to split pooled bams into individual
	bams dependent on output produced by vireo. </dd>
	<dt> demultiplex_no_argp.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules for demultiplexing using hashsolo and/or vireoSNP output and create a count matrix file. </dd>
	<dt> identify_swaps.snkmk </dt>
	<dd> This snakemake sub-workflow contains rules for identifying swaps using QTLtools-mbv. </dd>

<sup id="#ft1"><small>*</small></sup><small> - Consolidating into one</small>

<sup id="#ft2"><small>#</small></sup><small> - Not yet implemented</small>
</dl>
