# scRNA_snRNA_Pipeline
Snakemake Pipeline with "selectable" modules

## Requirements
This pipeline depends on the following packages/programs:
<ul>
	<li><a href="https://scanpy.readthedocs.io/en/stable/">Scanpy Manual</a></li>
	<li><a href="https://snakemake.readthedocs.io/en/stable/">Snakemake Manual</a></li>
	<li><a href="https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md">STARsolo Manual</a> </li>
	<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics">GC bias metrics Info</a></li>
	<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics">RNA seq metrics Info</a></li>
	<li><a href="https://github.com/pachterlab/kite">KITE <i>(kallisto indexing and tag extraction)</i></a></li>
	<li><a href="https://cellsnp-lite.readthedocs.io/en/latest/manual.html">cellSNP Manual</a></li>
	<li><a href="https://vireosnp.readthedocs.io/en/latest/manual.html">vireoSNP Manual</a></li>
	<li><a href="https://github.com/calico/solo#how-to-demultiplex-cell-hashing-data-using-hashsolo-cli">hashsolo Info</a></li>
</ul>

### Packages installed through conda
anaconda3/2018.12
Python 3.9.5

### Packages installed through pip
For R version 4.1.0

## Overview of the pipeline
Directed Acyclic Graph of the whole pipeline:
![DAG](images/Whole_pipeline.png)

## Available Modules

<ul>
<li> all</li>
<li> STARsolo</li>
<li> STARsolo_rnaseqmet</li>
<li> STARsolo_gcbiasmet</li>
<li> STARsolo_kb_solo</li>
<li> STARsolo_PICARD</li>
<li> STARsolo_gt_demux</li>
</ul>

### Module Specifications
<dl>
	<dt>Whole Pipeline</dt>
	<dd>This module includes alignment through STARsolo and PICARD's both programs (GcBiasMetrics and RnaSeqMetrics)</dd>
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


