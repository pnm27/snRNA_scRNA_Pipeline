#!/usr/bin/env python3

import os, glob2
from .utils import read_files_ext, get_demux_type
from .modality import get_modality
from .aligner import bam_paths, count_matrix as aligner_count_matrix
from .config_models import (
    demux_run, final_paths, final_sources, raw_output, swap_correction
)



# cellranger ----------------------------------------------------------------------------
def get_fastqs(wildcards, config):
    temp_fold=f"{config['fold_struct']}".format(**wildcards) # Add wildcards
    all_files = []
    all_files.extend(glob2.glob(f"{config['cDNA_fastqs_dir']}"
                f"{temp_fold}*{config['multiome_suffix']}"))
    all_files.extend(glob2.glob(f"{config['ATAC_fastqs_dir']}"
                f"{temp_fold}*{config['multiome_suffix']}"))

    return sorted(all_files)

# ---------------------------------------------------------------------------------------

# STARsolo ------------------------------------------------------------------------------
def get_r1_fastqs(wildcards, config):
    temp_fold=f"{config['fold_struct']}".format(**wildcards) # Add wildcards
    r1_files = sorted(glob2.glob(f"{config['cDNA_fastqs_dir']}{temp_fold}*{config['R1_suffix']}"))
    return r1_files


def get_r2_fastqs(wildcards, config):
    temp_fold=f"{config['fold_struct']}".format(**wildcards)
    r2_files = sorted(glob2.glob(f"{config['cDNA_fastqs_dir']}{temp_fold}*{config['R2_suffix']}"))
    return r2_files

# ---------------------------------------------------------------------------------------
# What a per-pool vcf (4-column 'vcf_info') is used for: the cellSNP
# pileup ('pileup') and/or vireo's donor genotypes ('vireo')
def _pool_vcf_use(config):
    use = config['gt_demux_pipeline'].get('pool_vcf_use', ['pileup', 'vireo'])
    return [u.lower() for u in ([use] if isinstance(use, str) else use or [])]

# Columns for a vcf_info file when provided with a vcf per pool
col_set = [
    "pool", 
    "n_dons", 
    "donors", # NEW ADDITION (comma-sep donor list)
    "vcf",
]
def previous_solo_h5ad(wildcards, config):
    """hashsolo h5ad (Rule A) of 'previous_run_tag', which add_vireo adds to."""
    return demux_run(
        config, 'solo', config['demultiplex'].get('previous_run_tag')
    ).h5ad


def get_cellsnp_inputs(wildcards, config):
    ret_dict={}
    ret_dict['barcodesFile'] = (
        f"{config['gt_demux_pipeline']['inp_for_cellsnp_dir']}"
        f"{config['fold_struct_filt_bc']}.txt"
    )
    ret_dict.update(bam_paths(config, get_modality(wildcards)))

    # 'cellsnp_ref_snps' is useful when per-Pool vcfs need to subset for SNVs
    if config['gt_demux_pipeline']['vcf_info'] is not None \
        and os.path.isfile(config['gt_demux_pipeline']['vcf_info']):
            samp_name = wildcards.pool # WILDCARDS
            temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
            n_cols = temp_df.shape[1]
            # For demux run wo_gt, we generally won't have per-donor vcf
            # Similarly, if there are just 2 columns and if it's not 
            # the demultiplexing run without genotype (like projects 
            # with only SNP array data) 
            if config['demux_run_type'] == 'wo_gt' or (
                n_cols == 2 and config['demux_run_type'] != 'wo_gt'
                ):
                ret_dict['regionsFile'] = config['gt_demux_pipeline']['genome_1k_ref']

            # If there are more than 3 columns then columns 4 and above
            # have per-donor vcfs 
            # (the 'above' is useful for multi-vcf setup) - TODO
            elif n_cols == 4 and 'pileup' in _pool_vcf_use(config):
                temp_df = temp_df.iloc[:, :4]
                temp_df.columns = col_set
                ret_dict['regionsFile'] = temp_df.loc[
                    temp_df["pool"] == samp_name, "vcf"].values[0]
            # Default
            else:
                ret_dict['regionsFile'] = config['gt_demux_pipeline']['genome_1k_ref']

    return ret_dict

# Consolidate this function and ret_dons
# MAKE THIS COMPATIBLE WITH NPSAD INPUTS WHERE WE HAVE MULTIPLE VCF POOL RUNS i.e. WITH FUNCTION
# MULTI_VCF
# i.e. An input to vireo rule (config['gt_demux_pipeline']['vcf_info']) should be either:
# 1) A three columned (no-header) file containing pool name, expected donors,
#    and corresponding vcf file (can be multiple columns) while config['gt_demux_pipeline']['donors_per_pool_file'] should be None
# 2) A 2-columned (no-header) txt or csv or tsv file with pool name and expected donors
#    and the path to a (single) vcf file (bgzf) containing all donors in config['gt_demux_pipeline']['donors_per_pool_file']
# 3) None, if the runs are to be without any vcf files and config['gt_demux_pipeline']['donors_per_pool_file'] should be 
#    a 2-columned (no-header) txt or csv or tsv file with pool name and expected donors
def get_vir_inputs(wildcards, config):
    ret_dict = {}
    # 'reuse_cellsnp_dir' points vireo to the cellSNP outputs of a previous
    # run (e.g. a wo_gt run), so cellSNP isn't re-run when SNPs haven't changed
    cellsnp_dir = (
        config['gt_demux_pipeline'].get('reuse_cellsnp_dir')
        or config['gt_demux_pipeline']['cellsnp_dir']
    )
    ret_dict['cellsnpCells'] = (
        f"{cellsnp_dir}"
        f"{config['fold_struct_gt_demux']}"
        f"{config['gt_demux_pipeline']['cellsnp_cells']}"    
    )

    if config['gt_demux_pipeline']['vcf_info'] is not None \
        and os.path.isfile(config['gt_demux_pipeline']['vcf_info']):
        temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
        n_cols = temp_df.shape[1]
        # Make snakemake's wildcard same as the value in the "pool" column
        samp_name = wildcards.pool # WILDCARDS
        # For condn. 2) - a single vcf containing all donors
        if n_cols == 2 and config['demux_run_type'] != 'wo_gt':
            ret_dict['donorFile'] = (
                config['gt_demux_pipeline']['donors_per_pool_file'] or []
            )
        # For condn. 1)
        elif n_cols == 4 and 'vireo' in _pool_vcf_use(config):
            temp_df = temp_df.iloc[:, :n_cols]
            temp_df.columns = col_set
            ret_dict['donorFile'] = (
                temp_df.loc[temp_df["pool"] == samp_name, "vcf"]
                .values[0]
                )
        else:
            ret_dict['donorFile'] = []  

    return ret_dict


# Demultiplex ---------------------------------------------------------------------------

def count_matrix(wildcards, config):
    """Count matrix (mtx) of the alignment, which Rule A starts from."""
    return aligner_count_matrix(config)


def get_inputs_demux_run(wildcards, config, tool):
    """Inputs of Rule A ('demux_run_<tool>'), see lib/config_models.py."""
    ret = {'mtx': count_matrix(wildcards, config)}

    if tool == 'solo':
        ret['calico_solo'] = raw_output(config, 'solo')

    elif tool == 'vireo':
        vireo_out = raw_output(config, 'vireo')
        if 'multiome' in config['last_step'].lower():
            # One vireo output per modality, in the order of the
            # '--suffix atac cdna' passed by 'get_params_demux'
            ret['vireo'] = [
                vireo_out.replace('{modality}', m) for m in ['ATAC', 'cDNA']
            ]
        else:
            ret['vireo'] = vireo_out
        if config['gt_demux_pipeline']['donorName_conv']['file'] is not None:
            ret['converter'] = config['gt_demux_pipeline']['donorName_conv']['file']

    return ret


def get_inputs_demux_final(wildcards, config):
    """Inputs of Rule B ('demux_final'): the Rule A outputs picked for the pool."""
    sources = final_sources(config, wildcards.pool) # WILDCARDS
    ret = {
        'h5ad': [s.h5ad for s in sources],
        'stats': [s.stats for s in sources if s.stats],
    }
    swap_file = swap_correction(config).get('file')
    if swap_file:
        ret['swap'] = swap_file
    return ret


# ---------------------------------------------------------------------------------------
# Split bam -----------------------------------------------------------------------------
def split_tool(config) -> str:
    """Tool of 'split_bams_pipeline: split_by: demux' (solo or vireo)."""
    demux = config['split_bams_pipeline']['split_by']['demux'].lower()
    if demux in ['vireo', 'vs']:
        return 'vireo'
    if demux in ["cs", "calico", "calico_solo", "hashsolo"]:
        return 'solo'
    raise ValueError(
        "Unknown 'split_bams_pipeline: split_by: demux' value: "
        f"'{demux}'"
    )


def get_bam(wildcards, config, suffix=None):
    """bam and bai of the pool's modality (suffix: see lib/modality.py)."""
    return bam_paths(config, get_modality(wildcards), suffix)


def bam_to_split(wildcards, config) -> str:
    """Which bam is split into per-donor bams.

    Genotype checks (gt_check) use the whole bam ('whole'), or the one
    subset to 'subset_chr' ('chr_subset'); otherwise the bam filtered by
    cell barcode ('cell_filtered').
    """
    if not config['gt_check']:
        return 'cell_filtered'
    if config['split_bams_pipeline']['subset_chr'] is None:
        return 'whole'
    return 'chr_subset'

# ---------------------------------------------------------------------------------------

# identify_swaps ------------------------------------------------------------------------

def get_bam_inputs(wildcards, config):
    
    if config['identify_swaps']['mbv_inp'] == 'vireo_outs':
        return (
            f"{config['split_bams_pipeline']['split_bams_dir']}"
            f"{config['fold_struct_bam_split2']}"
            f"{config['fold_struct_gt_demux_redo']}.bam"
        )
    else:
        return config['identify_swaps']['mbv_inp']

