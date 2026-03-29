#!/usr/bin/env python3

from .utils import read_files_ext


def get_params(wildcards, input, config, global_vars):
    params_dict = {
        "col_name": [config['gt_demux_pipeline']['demux_col'], "-c"],
        "bc_len": [config['gt_demux_pipeline']['barcode_len'], "-b"],
        "doub": [config['gt_demux_pipeline']['doublet'], "-d"],
        "neg": [config['gt_demux_pipeline']['negative'], "-n"],
        "na": [config['gt_demux_pipeline']['na'], "-e"],
        "max_mito": [config['max_mito_percentage'], "-m"],
        "min_genes": [config['min_genes_per_cell'], "-g"],
        "min_cells": [config['min_cells_per_gene'], "--min_cells"],
        "genes_info": [config['gene_info_file'], "--id2name"],
        "mito_prefix": [config['mito_prefix'], "--mito_prefix"],
 
    }
    ret_str = ''

    for k, v in params_dict.items():
        if k == "na" or k == "mito_prefix":
            ret_str += f'{v[1]} \'{v[0]}\' '
        else:
            ret_str += f'{v[1]} {v[0]} '

    if global_vars.ADD_SOLO or global_vars.ADD_VIREO or input[0].endswith('.h5ad'):
        ret_str += '--prev '
    
    if 'multiome' in config['last_step'].lower():
        ret_str += '--keep_barcode_suffix '
    
    if config['gt_demux_pipeline']['include_all_cells']:
        ret_str += '--keep_all_cells '

    return ret_str


def get_cmd_str_csnp(wildcards, input):
    ret_str = ""
    for x,y in zip(["-b", "-s"], input):
        ret_str+=f" {x} {y}"

    return ret_str



def get_cmd_str_vireo(wildcards, input, config) -> None | int:
    temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
    n_cols = 3 if config['gt_demux_pipeline']['vcf_info_columns']['vcf'] is not None else 2
    # Make snakemake's wildcard same as the value in the "pool" column
    samp_name = '-'.join(wildcards.pool.split('-')[:-1]) # WILDCARDS
    pool_col = config['gt_demux_pipeline']['vcf_info_columns']['pool']
    don_col = config['gt_demux_pipeline']['vcf_info_columns']['n_dons']
    ret_str = ""

    for x,y in zip(["-c", "-d"], input):
        ret_str+=f" {x} {y}"

    
    col_set1 = ["pool", "n_dons", "vcf"]
    col_set2 = ["pool", "n_dons"]
    # n_cols = ret_cols(config['gt_demux_pipeline']['vcf_info'])
    # if n_cols == 3:
    #     temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'], names=col_set1)
    # elif n_cols == 2:
    #     temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'], names=col_set2)
    # else:
    #     raise ValueError("Unexpected number of columns for 'vcf_info' file in 'gt_demux_pipeline'!!!")
    
    if samp_name.lower() in temp_df[pool_col].str.lower().values:
        ret_str += " -N " + str(temp_df.loc[temp_df[pool_col].str.lower() == samp_name.lower(), don_col].values[0])

    return ret_str
        

def get_umiTag(wildcards, config):
    # For ATAC return None
    if 'multiome' in config['last_step'].lower() and \
     'cdna' not in wildcards.pool.lower():
        return None
    else:
        config['gt_demux_pipeline']['umi_tag']