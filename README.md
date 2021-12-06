# scRNA_snRNA_Pipeline
Snakemake Pipeline with "selectable" modules

## Requirements
This pipeline executes the following set of modules:



## Selectable Modules
Directed Acyclic Graph of the whole pipeline:
![DAG](images/Whole_pipeline.png)

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
	<dd>This module includes alignment through STARsolo, PICARD's both programs (GcBiasMetrics and RnaSeqMetrics)</dd>
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

