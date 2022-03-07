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

#Parse wildcards from the file specified in config.yaml
round_num=[] # wildcard 'num'
sample_name=[] # wildcard 'id1'


# Limitting Step for the run of Snakemake, creating wildcards
with open(config['select_fastqs']) as fq:
    for line in fq:
        line_sp = line.split('/')
        round_num.append(line_sp[0])
        sample_name.append(line_sp[2].strip().replace('-cDNA', ''))



include: "rules/produce_targets.snkmk"
include: "rules/helper_functions.snkmk"
include: "rules/resources.snkmk"
include: "rules/calico_solo_demux.snkmk"
include: "rules/kite.snkmk"
# include: "rules/pheno_demux.snkmk"
include: "rules/picard_metrics.snkmk"
# include: "rules/split_fastqs.snkmk"
# include: "rules/split_fastqs_2.snkmk"
include: "rules/STARsolo.snkmk"
# include: "rules/test_split_bams.snkmk"
# For just gt purposes
include: "rules/split_bams_gt.snkmk"


rule all:
    input:
        produce_targets(conf_f=config)




