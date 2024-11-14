# snRNA_scRNA_Pipeline Introduction
## TODO:

- Miscellaneous:
   - [x] Add PICARD option in new_config file.
   - [ ] Write down schemas.
   - Add tutorials.
     - pooled snRNA seq
       - [ ] single wildcard
       - [ ] multiple wildcards
     - scRNA seq
       - [ ] single wildcard
       - [ ] multiple wildcards
     - [ ] Double HTOs
   - [ ] Remove dependency on STARsolo as an aligner.
   - [ ] For rules that use **genefull_matrices** make input function that take either *Gene* or *GeneFull* dependent on the project.
   - [ ] Combine sub-workflows split_bams and split_bams_gt.
      - [ ] Search Ranking of readthedocs (using config file for this too).
   - [ ] Might incorporate git submodules for repos on git that I use.
   - [ ] Add new Picard metrics.
   - [ ] Add options in config file to allow adding extra params for every software:
   - [ ] For reruns of vireo, provide a way to retain those information in update logs file.
   - [ ] Fix demultiplex_helper_funcs.py for double HTOs in function parse_file.
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
    - scripts:
- produce_targets:
    - new_config params:
    - snakemake_rules:
        - [ ] Simplify target functions.
    - scripts:
- resources:
    - scripts:
    - [ ] Add conditions for time and mem
- split_bams_gt*:
    - new_config params:
    - snakemake_rules:
        - [ ] Remove run directive
        - [ ] Input function that removes low mito cells.
    - scripts:
- split_bams*:
    - new_config params:
    - snakemake_rules:
        - [ ] Remove run directive
        - [ ] Input function that removes low mito cells.
    - scripts:
- STARsolo:
    - new_config params:
    - snakemake_rules:
        - [x] Remove run directive
        - [ ] WASP mode
    - scripts:
        - [ ] WASP mode
- Changelog:
    - Single wildcard is called now *pool* (Earlier mixed use of *num*, *id1* and *id2*).
    - Retained use of double wildcards.
    - Revised split_bams script:
        - [x] Consolidated gt and non-gt versions.
        - [x] Now, mito file is in params (earlier was an output)
        - [x] Now, bed file is in params (earlier was an inputs)
        - [x] Streamlined
    - Major revisions to **create_per_donor_bams.bash** script
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
    - Major revisions to **run_update_logs.sh** script
      - [x] Builds command from input using asssociative arrays.
      - [x] To emulate missingness (picard and/or demultiplexing) just use *empty* values.
      - [x] Now support for STARsolo 2.7.10 with **Final_out_MAP_2_7_10a_latest.tsv**.
    - Major revisions to **update_logs.py** script
      - [x] All optional parameters (except map_file, output_file, and bam_dir) expects one value or becomes *None* in its absence (with no argument value other values are used).
      - [x] To emulate missingness (picard and/or demultiplexing) just use *empty* values.
      - [x] Missingness of *picard_dir* implies not collecting GCBias and RNASeq Metrics.
      - [x] Similarly, missingness of *demul_dir* implies not demultiplexing info.
    - New file - **Final_out_MAP_2_7_10a_latest_info.xlsx** - contains more info related to **Final_out_MAP_2_7_10a_latest.tsv**.


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


## Changelog:

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
  -  Make the functions similar for demultiplexing with any method.
  - Fix issue with reading old wet_lab_info file to update (extension issues).
  - Some issue with create_wet_lab_info.py file (it misses to add some lines from certain files - try AMP ones)
  
  
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
