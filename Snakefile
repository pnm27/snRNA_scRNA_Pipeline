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


def get_all_inputs(conf_f, mm):
    if mm and (config['last_step'].endswith('.yaml') or config['last_step'].endswith('.yml')):
        # Validating config['last_step']
        # with open(config['last_step']) as fout:
        #     modules_df = pd.json_normalize(yaml.load(fout, Loader=yaml.SafeLoader))
        # validate(modules_df, "modules.schema.json")
        temp_l=[]
        with open(config['last_step']) as fout:
            modules_dict = yaml.load(fout, Loader=yaml.SafeLoader)

        # Load folder_structures yaml file
        with open(config['folder_structures']) as fout:
            fold_structs_dict = yaml.load(fout, Loader=yaml.SafeLoader)

        for k, v in modules_dict.items():
            temp_l.extend(produce_targets(conf_f=conf_f, last_step=v, wc_d=wildcards_list[k], fs=fold_structs_dict[k]))

        return temp_l

    elif not mm and not (config['last_step'].endswith('.yaml') or config['last_step'].endswith('.yml')):
        # Load folder_structures yaml file
        with open(config['folder_structures']) as fout:
            fold_structs_dict = yaml.load(fout, Loader=yaml.SafeLoader)

        return produce_targets(conf_f=conf_f, last_step=config['last_step'], wc_d=wildcards_list[k], fs=fold_structs_dict[k])

    # When a yaml file for modules is given but with only one module ( only for practical use when checking a part of multi_module setup)
    elif not mm and (config['last_step'].endswith('.yaml') or config['last_step'].endswith('.yml')):
        with open(config['last_step']) as fout:
            modules_dict = yaml.load(fout, Loader=yaml.SafeLoader)

        # Load folder_structures yaml file
        with open(config['folder_structures']) as fout:
            fold_structs_dict = yaml.load(fout, Loader=yaml.SafeLoader)

        for k,v in modules_dict.items():
            return produce_targets(conf_f=conf_f, last_step=v, wc_d=wildcards_list[k], fs=fold_structs_dict[k])



wildcard_constraints:
    f_s: "[A-Za-z]+(?=/)?"



rule all:
    input:
        get_all_inputs(conf_f=config, mm=MULTI_MODULES)
        




