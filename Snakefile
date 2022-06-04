import os
#from collections import OrderedDict
import glob2, re, math
import pandas as pd
from helper_py_scripts import update_logs
from snakemake.utils import validate
from itertools import repeat
from snakemake.utils import min_version


min_version("6.4.0")


configfile: "new_config.yaml"
#validate(config, "config.schema.json")  # Path to the scefic schema


include: "rules/input_processing.snkmk" # Process input files to create lists for wildcards
include: "rules/produce_targets.snkmk" # Produce target files for the pipeline using previously created lists for wildcards for the input files
include: "rules/helper_functions.snkmk"
include: "rules/resources.snkmk"
include: "rules/STARsolo.snkmk"
include: "rules/picard_metrics.snkmk"
include: "rules/kite.snkmk"
include: "rules/calico_solo_demux.snkmk"
# include: "rules/pheno_demux.snkmk"
# include: "rules/split_fastqs.snkmk"
# include: "rules/split_fastqs_2.snkmk"
include: "rules/split_bams.snkmk"
# For just gt purposes
# include: "rules/split_bams_gt.snkmk"
# include: "rules/pheno_demux2.snkmk"


rule all:
    input:
        produce_targets(conf_f=config)




