"""
Author: Prashant N M
Affiliation: Mount Sinai School of Medicine, Department of Psychiatry
Aim: Snakemake workflow for scRNAseq and snRNAseq supporting multiplexed pools
Date: <mod_date>
Run: indirectly run through run_snakemake.sh
Latest modification: 
  - Adding multi_module support
Further Plans:
  - generalize aligner
  - add 'message' for each rule
  - add 'log' for each rule
  - add 'report' to the pipeline
"""



import os
import glob2, re, math
import pandas as pd
from snakemake.utils import validate
from itertools import repeat
from snakemake.utils import min_version


min_version("6.4.0")


configfile: "new_config.yaml"
#validate(config, "config.schema.json")  # Path to the specific schema

# REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Set few global variables (DON'T CHANGE FOR INITIAL RUNS (FOR EACH PROJECT)!)
# FOr subsequent runs on the same project i.e. when there exist previous
# outputs of either demultiplex runs: CHECK THE DOCUMENTATION FOR HOW TO 
# CHANGE THESE VARS

class global_vars():
  ONLY_SOLO = False
  ONLY_VIREO = False
  BOTH_DEMUX = False # Not yet implemented the rule
  ADD_SOLO = False # When a demultiplex run with vireoSNP has been done
  ADD_VIREO = False # When a demultiplex run with calico_solo has been done
  SPLIT_BY_SOLO = False # When subset pooled bams by output of calico_solo
  SPLIT_BY_VIREO = True # When subset pooled bams by output of vireo


include: "rules/input_processing.smk"
include: "rules/produce_targets.smk"
include: "rules/helper_functions.smk"
include: "rules/resources.smk"
include: "rules/STARsolo.smk"
include: "rules/picard_metrics.smk"
include: "rules/kite.smk"
include: "rules/calico_solo_demux.smk"
# include: "rules/pheno_demux3.smk" # DEPRACATED NAME
include: "rules/genotype_demux.smk"
include: "rules/demultiplex.smk"
include: "rules/split_bams.smk"
include: "rules/identify_swaps.smk"
include: "rules/deeptools.smk"
include: "rules/cellranger.smk"


rule all:
    input:
        produce_targets(conf_f=config, last_step=config['last_step'], wc_d=wildcards_list)
        




