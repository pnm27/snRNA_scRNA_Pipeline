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
#from collections import OrderedDict
import glob2, re, math
import pandas as pd
from snakemake.utils import validate
from itertools import repeat
from snakemake.utils import min_version


min_version("6.4.0")


configfile: "new_config.yaml"
#validate(config, "config.schema.json")  # Path to the scefic schema

# Set few global variables (DON'T CHANGE!)
ONLY_SOLO = False
ONLY_VIREO = False
BOTH_DEMUX = False # Not eyt implemented the rule
ADD_SOLO = False # When a demultiplex run with vireoSNP has been done
ADD_VIREO = False # When a demultiplex run with calico_solo has been done, modify it later


include: "rules/input_processing.snkmk"
include: "rules/produce_targets.snkmk"
include: "rules/helper_functions.snkmk"
include: "rules/resources.snkmk"
include: "rules/STARsolo.snkmk"
include: "rules/picard_metrics.snkmk"
include: "rules/kite.snkmk"
include: "rules/calico_solo_demux.snkmk"
# include: "rules/pheno_demux.snkmk" # GOD knows why
# include: "rules/split_fastqs.snkmk" # outdated
# include: "rules/split_fastqs_2.snkmk" # outdated
include: "rules/pheno_demux3.snkmk"
include: "rules/demultiplex_no_argp.snkmk"
include: "rules/split_bams.snkmk"
# For just gt purposes
# include: "rules/split_bams_gt.snkmk"
# include: "rules/pheno_demux2.snkmk"


rule all:
    input:
        produce_targets(conf_f=config)
        




