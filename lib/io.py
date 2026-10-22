#!/usr/bin/env python3

import os, glob2
from .utils import read_files_ext, get_demux_type, is_genotyped
from .modality import get_modality, is_multiome
from .aligner import bam_paths, count_matrix as aligner_count_matrix
from .config_models import (
    demux_run, final_paths, final_sources, raw_output, swap_correction
)



# cellranger ----------------------------------------------------------------------------
def get_fastqs(wildcards, config):
    temp_fold=f"{config['modules']['STARsolo']['fold_struct']}".format(**wildcards) # Add wildcards
    all_files = []
    all_files.extend(glob2.glob(f"{config['experiment']['inputs']['fastqs_dir']['cDNA']}"
                f"{temp_fold}*{config['experiment']['inputs']['multiome_suffix']}"))
    all_files.extend(glob2.glob(f"{config['experiment']['inputs']['fastqs_dir']['ATAC']}"
                f"{temp_fold}*{config['experiment']['inputs']['multiome_suffix']}"))

    return sorted(all_files)

# ---------------------------------------------------------------------------------------

# STARsolo ------------------------------------------------------------------------------
def pool_vcf(wildcards, config):
    """Per-pool vcf of the pool (the 'vcf' column of a 4-column 'vcf_info')."""
    vcf_info = config['modules']['gt_demux']['vcf_info']
    if vcf_info is None or not os.path.isfile(vcf_info):
        raise ValueError(
            "STARsolo's WASP mode needs a per-pool vcf: set "
            "'modules: gt_demux: vcf_info' to a file with a vcf column "
            f"(currently: {vcf_info})"
        )
    temp_df = read_files_ext(vcf_info)
    if temp_df.shape[1] < 4:
        raise ValueError(
            f"STARsolo's WASP mode needs a vcf column in {vcf_info} "
            "(pool, n_donors, donors, vcf); it has no per-pool vcf"
        )
    temp_df = temp_df.iloc[:, :4]
    temp_df.columns = col_set
    vcf = temp_df.loc[temp_df["pool"] == wildcards.pool, "vcf"].dropna() # WILDCARDS
    if vcf.empty:
        raise ValueError(
            f"No vcf for the pool '{wildcards.pool}' in {vcf_info}, which "
            "STARsolo's WASP mode needs"
        )
    return vcf.values[0]


def wasp_vcf(config) -> str:
    """Pool's vcf for STAR's WASP mode (see rule 'wasp_vcf')."""
    c = config['modules']['STARsolo']
    return f"{c['bams_dir']}{config['modules']['STARsolo']['fold_struct']}_wasp.vcf"


def get_r1_fastqs(wildcards, config):
    temp_fold=f"{config['modules']['STARsolo']['fold_struct']}".format(**wildcards) # Add wildcards
    r1_files = sorted(glob2.glob(f"{config['experiment']['inputs']['fastqs_dir']['cDNA']}{temp_fold}*{config['experiment']['inputs']['R1_suffix']}"))
    return r1_files


def get_r2_fastqs(wildcards, config):
    temp_fold=f"{config['modules']['STARsolo']['fold_struct']}".format(**wildcards)
    r2_files = sorted(glob2.glob(f"{config['experiment']['inputs']['fastqs_dir']['cDNA']}{temp_fold}*{config['experiment']['inputs']['R2_suffix']}"))
    return r2_files

# ---------------------------------------------------------------------------------------
# What a per-pool vcf (4-column 'vcf_info') is used for: the cellSNP
# pileup ('pileup') and/or vireo's donor genotypes ('vireo')
def _pool_vcf_use(config):
    use = config['modules']['gt_demux'].get('pool_vcf_use', ['pileup', 'vireo'])
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
        config, 'solo', config['modules']['demultiplex'].get('previous_run_tag')
    ).h5ad


def get_cellsnp_inputs(wildcards, config):
    ret_dict={}
    ret_dict['barcodesFile'] = (
        f"{config['modules']['gt_demux']['inp_for_cellsnp_dir']}"
        f"{config['modules']['gt_demux']['fold_struct_filt_bc']}.txt"
    )
    ret_dict.update(bam_paths(config, get_modality(wildcards)))

    # 'cellsnp_ref_snps' is useful when per-Pool vcfs need to subset for SNVs
    if config['modules']['gt_demux']['vcf_info'] is not None \
        and os.path.isfile(config['modules']['gt_demux']['vcf_info']):
            samp_name = wildcards.pool # WILDCARDS
            temp_df = read_files_ext(config['modules']['gt_demux']['vcf_info'])
            n_cols = temp_df.shape[1]
            # For demux run wo_gt, we generally won't have per-donor vcf
            # Similarly, if there are just 2 columns and if it's not 
            # the demultiplexing run without genotype (like projects 
            # with only SNP array data) 
            if not is_genotyped(config) or n_cols == 2:
                ret_dict['regionsFile'] = config['modules']['gt_demux']['genome_1k_ref']

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
                ret_dict['regionsFile'] = config['modules']['gt_demux']['genome_1k_ref']

    return ret_dict

# Consolidate this function and ret_dons
# MAKE THIS COMPATIBLE WITH NPSAD INPUTS WHERE WE HAVE MULTIPLE VCF POOL RUNS i.e. WITH FUNCTION
# MULTI_VCF
# i.e. An input to vireo rule (config['modules']['gt_demux']['vcf_info']) should be either:
# 1) A three columned (no-header) file containing pool name, expected donors,
#    and corresponding vcf file (can be multiple columns) while config['modules']['gt_demux']['donors_per_pool_file'] should be None
# 2) A 2-columned (no-header) txt or csv or tsv file with pool name and expected donors
#    and the path to a (single) vcf file (bgzf) containing all donors in config['modules']['gt_demux']['donors_per_pool_file']
# 3) None, if the runs are to be without any vcf files and config['modules']['gt_demux']['donors_per_pool_file'] should be 
#    a 2-columned (no-header) txt or csv or tsv file with pool name and expected donors
def get_vir_inputs(wildcards, config):
    ret_dict = {}
    # 'reuse_cellsnp_dir' points vireo to the cellSNP outputs of a previous
    # run (e.g. a wo_gt run), so cellSNP isn't re-run when SNPs haven't changed
    cellsnp_dir = (
        config['modules']['gt_demux'].get('reuse_cellsnp_dir')
        or config['modules']['gt_demux']['cellsnp_dir']
    )
    ret_dict['cellsnpCells'] = (
        f"{cellsnp_dir}"
        f"{config['modules']['gt_demux']['fold_struct']}"
        f"{config['modules']['gt_demux']['cellsnp_cells']}"    
    )

    if config['modules']['gt_demux']['vcf_info'] is not None \
        and os.path.isfile(config['modules']['gt_demux']['vcf_info']):
        temp_df = read_files_ext(config['modules']['gt_demux']['vcf_info'])
        n_cols = temp_df.shape[1]
        # Make snakemake's wildcard same as the value in the "pool" column
        samp_name = wildcards.pool # WILDCARDS
        # For condn. 2) - a single vcf containing all donors
        if n_cols == 2 and is_genotyped(config):
            ret_dict['donorFile'] = (
                config['modules']['gt_demux']['donors_per_pool_file'] or []
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
        if is_multiome(config):
            # One vireo output per modality, in the order of the
            # '--suffix atac cdna' passed by 'get_params_demux'
            ret['vireo'] = [
                vireo_out.replace('{modality}', m) for m in ['ATAC', 'cDNA']
            ]
        else:
            ret['vireo'] = vireo_out
        if config['modules']['gt_demux']['donorName_conv']['file'] is not None:
            ret['converter'] = config['modules']['gt_demux']['donorName_conv']['file']

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
    """Tool of 'modules: split_bams: split_by: demux' (solo or vireo)."""
    demux = config['modules']['split_bams']['split_by']['demux'].lower()
    if demux in ['vireo', 'vs']:
        return 'vireo'
    if demux in ["cs", "calico", "calico_solo", "hashsolo"]:
        return 'solo'
    raise ValueError(
        "Unknown 'modules: split_bams: split_by: demux' value: "
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
    if not config['workflow']['gt_check']:
        return 'cell_filtered'
    if config['modules']['split_bams']['subset_chr'] is None:
        return 'whole'
    return 'chr_subset'

# ---------------------------------------------------------------------------------------

# identify_swaps ------------------------------------------------------------------------

def get_bam_inputs(wildcards, config):
    
    if config['modules']['identify_swaps']['mbv_inp'] == 'vireo_outs':
        return (
            f"{config['modules']['split_bams']['split_bams_dir']}"
            f"{config['modules']['split_bams']['fold_struct_bams']}"
            f"{config['modules']['split_bams']['fold_struct_donor']}.bam"
        )
    else:
        return config['modules']['identify_swaps']['mbv_inp']

