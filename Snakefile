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



import sys
from snakemake.utils import validate
from snakemake.utils import min_version


min_version("8.25.0")

assert sys.version_info >= (3, 10), "This pipeline needs python version >= 3.10!"


configfile: "new_config.yaml"
#validate(config, "config.schema.json")  # Path to the specific schema


# DEPRACATED - USE new_config.yaml 'demux_type'
# class global_vars():
#   ONLY_SOLO = False
#   ONLY_VIREO = False
#   BOTH_DEMUX = False 
#   ADD_SOLO = False # When a demultiplex run with vireoSNP has been done
#   ADD_VIREO = False # When a demultiplex run with calico_solo has been done
#   SPLIT_BY_SOLO = False # When subset pooled bams by output of calico_solo
#   SPLIT_BY_VIREO = True # When subset pooled bams by output of vireo


include: "rules/input_processing.smk"
include: "rules/produce_targets.smk"
include: "rules/resources.smk"
include: "rules/STARsolo.smk"
include: "rules/picard_metrics.smk"
include: "rules/kite.smk"
include: "rules/calico_solo_demux.smk"
include: "rules/genotype_demux.smk"
include: "rules/demultiplex.smk"
include: "rules/split_bams.smk"
include: "rules/identify_swaps.smk"
include: "rules/deeptools.smk"
include: "rules/cellranger.smk"


rule all:
    input:
        produce_targets


# dummy rule to accommodate checkpoints output
# Gets triggered by presence of any of 'split_bams' or 
# 'identify_swaps' in last_step
rule resolve_pool_targets:
    input:
        produce_targets_dynamic
    output:
        "resolved/{pool}/{modality}.done"
    shell:
        "touch {output}"
        




