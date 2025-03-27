# snRNA_scRNA_Pipeline Introduction

## TODO

- Miscellaneous:
  - [ ] Write down schemas.
  - [ ] Add support for multiome in **update_logs.py**.
  - Add tutorials.
    - pooled snRNA seq
      - [ ] single wildcard
      - [ ] multiple wildcards
    - scRNA seq
      - [ ] single wildcard
      - [ ] multiple wildcards
    - Double HTOs
  - Generalize Outputs for split bam pipeline (whole of **split_bams.smk**).
  - [ ] Remove dependency on STARsolo as an aligner.
  - [ ] For rules that use **genefull_matrices** make input function that take either *Gene* or *GeneFull* dependent on the project.
  - [x] Combine sub-workflows split_bams and split_bams_gt.
    - [ ] Search Ranking of readthedocs (using config file for this too).
  - [ ] Might incorporate git submodules for repos on git that I use.
  - [ ] Add new Picard metrics.
  - [x] Add options in config file to allow adding extra params for every software:
  - [ ] For reruns of vireo, provide a way to retain those information in update logs file.
  - [ ] Retain information from log runs:
    - [ ] Use command like following to extract *variants matched to genotype VCF out of total variants piledup for the pool*.
  
      ```shell
      find -mindepth 1 -maxdepth 1 -type d -exec sh -c 'a=$(sed "s#\./##g" <<< {}); b=$(ls -ltr ${a} | tail -1 | rev | cut -d " " -f1 | rev); grep "variants matched to donor VCF" ${a}/${b}' \;
      ```

  - [ ] Fix demultiplex_helper_funcs.py for double HTOs in function parse_file.
  - [ ] Simplify structure of wildcards.
    - Folder structure should include
- analyse_vireo:
  - new_config params:
  - snakemake_rules:
  - scripts:
- calico_solo_demux:
  - new_config params:
  - snakemake_rules:
  - scripts:
- demultiplex:
  - new_config params:
  - snakemake_rules:
  - [ ] Employ a strategy for final count matrix dir (file dir in config file) for the cases:
    - when both demultiplex software are run simultaneously.
    - When there's an order (try to name each run separately or at least keep the order somewhere mentioned).
  - scripts:
    - [ ] when adding calico_solo or vireo include the demultiplex file (file containing demux stats) as input and append to it.
    - [ ] For reruns of vireo, provide a way to retain those information in demultiplex info file.
- helper_functions:
  - new_config params:
  - snakemake_rules:
  - scripts:
- identify_swaps:
  - new_config params:
  - snakemake_rules:
  - scripts:
- input_processing:
  - new_config params:
  - snakemake_rules:
  - scripts:
- kite:
  - new_config params:
  - snakemake_rules:
    - [x] Remove run directive
  - scripts:
- pheno_demux3:
  - new_config params:
  - snakemake_rules:
    - [x] Remove run directive
    - [ ] Beautify the function get_filt_barcodes.
  - scripts:
- picard_metrics:
  - new_config params:
  - snakemake_rules:
    - [ ] Add CollectInsertSizeMetrics for ATAC part of multiome.
  - scripts:
- produce_targets:
  - new_config params:
  - snakemake_rules:
    - [ ] Simplify target functions.
  - scripts:
- STARsolo:
  - new_config params:
  - snakemake_rules:
    - [x] Remove run directive
    - [ ] WASP mode
  - [ ] Issues with using wildcard **vcf_type** in the rule *demux_samples*.
  - [ ] Issues with output dir selection in the *demux_samples* i.e. automatically pick output dir.
  - [ ] Revamp wildcards so that varying output dirs are corrected accordingly:
    - [ ] Demux output i.e. if only one demux method needs to be used or simultaneously both.
    - [ ] splitting bams is for finalizing or genotype purposes.
- cellranger:
  - [ ] Support for cellranger.
    - Support for cellranger-arc count.
      - [x] Add support for alignment.
      - [ ] Add support for ATAC-based vireo demultiplexing.

This pipeline intends to not only make complex {term}`preprocessing` workflows easy (e.g. snRNA seq with pooled samples, double HTOs, etc.) but also to facilitate the use of common workflows used for preprocessing by providing *readymade* different combinations of softwares/tools (see {ref}`selectable <selectable-modules>` modules for more options). 

It also supports various software/pipeline for scRNA seq pre-processing.

The highlights of the pipeline are:
<ul>
  <li> Streamlined processes to modify parameters for each program through a single yaml file </li>
  <li> Easily modifiable to accomodate more rules </li>
  <li> Can be used for both individual samples as well as multiplexed pools </li>
  <li> Preserve folder structures (mirroring fastqs' folder structures) </li>
  <li> Organize outputs from each module </li>
  <li> Select multiple pre-set modules that simplifies usage across multiple projects </li>
</ul>

## Changelog

- Changed param name in demultiplex info from *Unique genes* to *gene_ids with an associated gene_name*.
- Added new param in demultiplex info file to add more stats when remove gene IDs without an associated gene name.
- Added an option to run cellSNP without any ref vcfs (1000 Genomes Project vcf is min requirement)
- Now create_wet_lab_info scripts can:
  - Run without a converter file
  - Save donor file along with the wet lab compilation file
  - argparse documented
- Fixed an issue with create_wet_lab_info.py file
- create_wet_lab_info.py file now mirrors actions for donor and multiplex compilations.
- Changed name of the rule demux_samples_calico_solo_STARsolo to demux_samples.
- Changed the *demux_info* parameter to optional (from positional) in demultiplex_no_argp.snkmk's rule that handles adding new demux to a final count matrix.
- Added working argparse to demul_samples_no_argp.py script.
- Changed the name of sub-workflow demultiplex_no_argp.snkmk to demultiplex.snkmk
- Changed the name of sub-workflow demul_samples_no_argp.py to demul_samples.py
- Fix demultiplex_no_argp.snkmk's rule that handles adding new demux to a final count matrix.
- Add an option (in config file) to create h5ads when demultiplexing (demultiplex_no_argp.snkmk) or not (can be used as switch when doing gt checks and finalizing donor assignment).
- Add an option for the rule cellSNP when ref SNPs vcf need not be subsetted further.
- Make the functions similar for demultiplexing with any method.
- Fix issue with reading old wet_lab_info file to update (extension issues).
- Some issue with create_wet_lab_info.py file (it misses to add some lines from certain files - try AMP ones)
- Single wildcard is called now *pool* (Earlier mixed use of *num*, *id1* and *id2*).
- Retained use of double wildcards.
- Revised split_bams script:
  - [x] Consolidated gt and non-gt versions.
  - [x] Now, mito file is in params (earlier was an output)
  - [x] Now, bed file is in params (earlier was an inputs)
  - [x] Streamlined
- Major revisions to *create_per_donor_bams.bash* script
  - [x] Consolidated gt and non-gt versions.
  - [x] Handles saving mito_file much elegantly.
  - [x] Supports argument parsing (with support for older positional args)
  - [x] Doesn't expect directories, provided as inputs, to follow logic - dirs should end with '/.
- Changed name of workflow from **pheno_demux3.snkmk** to **genotype_demux.snkmk**
- The rule create_inp_splitBams now:
  - [x] Uses a consolidated script to create the barcode files using both, h5ad and raw files.
  - [x] Removed the option to overwrite the outputs as before (but still present in python script).
  - [x] Rule supports single demux (while the script can handle multiple).
  - [x] h5ad input alone support for calico_solo while vireo output is supported as is.
- Major revisions to *run_update_logs.sh* script
  - [x] Builds command from input using asssociative arrays.
  - [x] To emulate missingness (picard and/or demultiplexing) just use *empty* values.
  - [x] Now support for STARsolo 2.7.10 with *Final_out_MAP_2_7_10a_latest.tsv*.
- Major revisions to *update_logs.py* script
  - [x] All optional parameters (except map_file, output_file, and bam_dir) expects one value or becomes *None* in its absence (with no argument value other values are used).
  - [x] To emulate missingness (picard and/or demultiplexing) just use *empty* values.
  - [x] Missingness of *picard_dir* implies not collecting GCBias and RNASeq Metrics.
  - [x] Similarly, missingness of *demul_dir* implies not demultiplexing info.
- New file - *Final_out_MAP_2_7_10a_latest_info.xlsx* - contains more info related to *Final_out_MAP_2_7_10a_latest.tsv*.
- Changed the section name from *demux_pipeline* to *hashsolo_demux_pipeline* in **new_config.yaml** file.
- **calico_solo_demux.smk**: At lines numbered 3 and 9.
- **create_logs.smk**: At lines numbered 11 and 52.
- **demultiplex.smk**: At lines numbered 10,12,32,34,72,74,78,80,124,125,132,135,138,328,331, and 333.
- **genotype_demux.smk**: At line numbered 45.
- **kite.smk**: At lines numbered 417 and 437.
- **produce_targets.smk**: At lines numbered 159-162.
- **split_bams.smk**: At lines numbered 10,16, and 18.
- Consolidated the section name from *split_bams_pipeline_gt_demux* in *split_bams_pipeline_gt* (earlier used for calico_solo based split bams) in **new_config.yaml** file.
  - **split_bams.smk**: At lines numbered 44,45,47 and 52.
  - **produce_targets.smk**: At line numbered 111.
  - **identify_swaps.smk**: At line numbered 3.
- Added support for Snakemake transition to version > 8.
  - Added *workflow_profile/config.yaml* which gets reflected in *run_snakemake.sh*.
  - To emulate previous behavior's for profile manually edited the *lsf_executor_plugin* and added ENV variable.
  - *lsf.yaml* still is present for snakemake \< v8.
  - Changed the *threads* directive and replaced with resources: *cpus_per_task*.
- Removed dependence on **resources.smk**. Instead all resource requirements are within each snakemake file.
- Simplified the rule *cellSNP* in **genotype_demux.smk**.
  - Now a parameter *cmd_str_csnp* function to replace the use of indexed arrays in shell (was working Snakemake \< 8).
  - Simplified the commandline for execution.
  - Now the rule picks the biggest of given *n_proc* in *new_config.yaml* and twice of the number of cpus provided.
- Simplified the rule *vireoSNP* in **genotype_demux.smk**.
  - Now a parameter *cmd_str_vireo* function to replace the use of indexed arrays in shell (was working Snakemake \< 8).
  - Now file specified in *vcf_info* (relates to the rule *vireoSNP*) in *new_config.yaml* is expected to have headers.
  - Simplified the commandline for execution.
- Usage of *pd.concat* now in concordance with [FutureWarning](https://github.com/pandas-dev/pandas/blob/a0babcb2c63dd721ea47e75f6229c5fe727b2395/pandas/core/internals/concat.py#L492) in *update_logs.py*.
- In **new_config.yaml**, changed 
  - *gt_conv* to *file* in *donorName_conv* in *gt_demux_pipeline*.
  - *mito* to *mito_prefix*. Reflected in **demultiplex.smk**, **split_bams.smk** and **calico_solo_demux.smk**
  - *gt_check* in *gt_check* to *gt_check*. Reflected in **split_bams.smk** and **produce_targets.smk**.
  - Added *demultiplex* section for the rule *demux_samples_both*.
- Removed mode='w+' when creating outputs in *create_Feat_Barc.py*.
- Added *multiome_alignment* as a new module. Created **cellranger.smk**, which currently support cellranger arc count only.
- Added multiome demultiplexing support for the following rules:
  - genotype_demux
    - Fixed UMItag selection in *cellSNP* for multiome-ATAC.
    - Fixed issues for multiome in *create_inp_cellSNP* (to add -1 in cell barcodes as bam by cellranger has "-1" suffix)
  - demultiplex
  - Change 'vcf_type' wildcard to support both multi-vcf and multiome setup.
- Added multiome support for splitting bams (using both cDNA and ATAC modalities).
- Simplified shell script in the rules *demux_samples* and *add_obs_to_final_count_matrix* in **demultiplex.smk**
- Added PICARD option in new_config file.
- In the module **demultiplex.smk**, changed logic of the rules *demux_samples* and *add_obs_to_final_count_matrix*:
  - 3 new rules reflect 3 different output types i.e. demux_samples_solo, demux_samples_vireo, and demux_samples_both.
  - Reflected in **produce_targets.smk**.
  - Now, support for multiome through the rule **demux_samples_vireo** (set global *ONLY_VIREO*).
  - Added support for multiome in *demul_samples.py* through append mode for vireo.
  
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
