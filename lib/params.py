#!/usr/bin/env python3

from .utils import read_files_ext


# Columns for a vcf_info file when provided with a vcf per pool
col_set = [
    "pool", 
    "n_dons", 
    "donors", # NEW ADDITION (comma-sep donor list)
    "vcf",
]

def get_params_cics(wildcards, config, input):
    params_dict = {
        "col_name": [config['gt_demux_pipeline']['demux_col'], "-c"],
        "bc_len": [config['gt_demux_pipeline']['barcode_len'], "-b"],
        "doub": [config['gt_demux_pipeline']['doublet'], "-d"],
        "neg": [config['gt_demux_pipeline']['negative'], "-n"],
        "na": [config['gt_demux_pipeline']['na'], "-e"],
        "max_mito": [config['max_mito_percentage'], "-m"],
        "min_genes": [config['min_genes_per_cell'], "-g"],
        "min_cells": [config['min_cells_per_gene'], "--min_cells"],
        "genes_info": [
            config['STARsolo_pipeline']['genome_pick']['gene_info_file'], 
            "--id2name"
            ],
        "mito_prefix": [config['mito_prefix'], "--mito_prefix"],
 
    }
    ret_str = ''

    for k, v in params_dict.items():
        if k == "na" or k == "mito_prefix":
            ret_str += f'{v[1]} \'{v[0]}\' '
        else:
            ret_str += f'{v[1]} {v[0]} '

    if config['demux_type'] in ['add_solo', 'add_vireo'] or input[0].endswith('.h5ad'):
        ret_str += '--prev '
    
    if 'multiome' in config['last_step'].lower():
        ret_str += '--keep_barcode_suffix '
    
    if config['gt_demux_pipeline']['include_all_cells']:
        ret_str += '--keep_all_cells '

    return ret_str


def get_procs_csnp(wildcards, config, resources):
    user_def = config['gt_demux_pipeline']['n_proc']
    expected_thresh = resources.cpus_per_task * 2
    return user_def if user_def > expected_thresh else expected_thresh


# DEPRACATED
# def get_cmd_str_csnp(wildcards, input):
#     ret_str = ""
#     inp_list = [
#         input.barcodefile,
#         input.bam
#     ]
#     for x,y in zip(["-b", "-s"], inp_list):
#         ret_str+=f" {x} {y}"

#     return ret_str



def get_cmd_str_vireo(wildcards, input, config) -> None | int:
    temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
    # Make snakemake's wildcard same as the value in the "pool" column
    samp_name = wildcards.pool # WILDCARDS
    ret_str = ""
    # for x,y in zip(["-c", "-d"], input):
    #     ret_str+=f" {x} {y}"
    if input['donorFile']:
        ret_str += f" -d {input['donorFile']}"

    temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
    cols = col_set[:2]
    temp_df = temp_df.iloc[:, :2]
    temp_df.columns = cols

    ret_str += " -N " + str(
        temp_df.loc[temp_df['pool'].str.lower() == samp_name.lower(), 
            'n_dons'].values[0]
        )

    return ret_str
        

def get_umiTag(wildcards, config):
    # For ATAC return None
    if 'multiome' in config['last_step'].lower() and \
     'cdna' not in wildcards.pool.lower():
        return None
    else:
        return config['gt_demux_pipeline']['umi_tag']


# 
def get_params_demux(wildcards, input, output, config):
    params_dict = {
        "demux_info": [output[1], "--demux_info"],
        "cols": [config['hashsolo_demux_pipeline']['columns_to_pick'], "--columns"],
        "wet_lab_file": [config['wet_lab_info'], "--wet_lab_file"],
        "hto_sep": [config['hashsolo_demux_pipeline']['hto_sep'], "--hto_sep"],
        "max_mito": [config['max_mito_percentage'], "-m"],
        "min_genes": [config['min_genes_per_cell'], "-g"],
        "min_cells": [config['min_cells_per_gene'], "--min_cells"],
        "mito_prefix": [config['mito_prefix'], "--mito_prefix"],
        "cs_conv": [config['hashsolo_demux_pipeline']['SubID_convert'], "--no-subid_convert"],
        "conv_f_nheaders": [
            config['gt_demux_pipeline']['donorName_conv']['header_lev'], 
            "--converter_file_headerNlev"],
        "conv_f_pool_col": [
            config['gt_demux_pipeline']['donorName_conv']['pool_col'], 
            "--conv_file_pool_column"],
        "conv_f_donor_col": [
            config['gt_demux_pipeline']['donorName_conv']['donor_col'], 
            "--conv_file_donor_column"],
        "conv_f_conv_col": [
            config['gt_demux_pipeline']['donorName_conv']['convert_col'], 
            "--conv_file_conv_column"],
        "new_h5ad_col": [
            config['gt_demux_pipeline']['donorName_conv']['new_h5ad_colname'], 
            "--h5ad_new_classify_colname"],
        "swap_correction_df": [
            config['gt_demux_pipeline']['swap_correction_df'],
            "--swap_correct",
        ]
    }
    gene_info_file = config['STARsolo_pipeline']['genome_pick']['gene_info_file']
    solo_inp = ["--calico_solo"]
    vireo_inp = ["--vireo_out", "--converter_file"]
    multiome_inp = ["--vireo_out", "--vireo_out", "--converter_file"]
    multiome_suffix = " --suffix atac cdna "
    inp_files = [""] # Keep first value empty as input[0] is positional arg
    multi_module = ['multiome', 'gt_demux']
    multi_module = all([ m in config['last_step'].lower() for m in multi_module])
    ret_str = ''
    pos_args = ''

    for k, v in params_dict.items():
        if k == "hto_sep" or k == "mito_prefix":
            ret_str += f'{v[1]} \'{v[0]}\' '
        elif k == 'cs_conv' and not v[0]:
            ret_str += f'{v[1]} '
        elif k == 'cs_conv' and v[0]:
            continue
        else:
            ret_str += f'{v[1]} {v[0]} '

    if config['demux_type'].lower() in ['solo', 'add_solo']:
        inp_files += solo_inp
    elif config['demux_type'].lower() in ['vireo', 'add_vireo']:
        if multi_module:
            inp_files += multiome_inp
        else:
            inp_files += vireo_inp
    elif config['demux_type'].lower() == 'both':
        inp_files += solo_inp + vireo_inp

    for k, v in zip(inp_files, input):
        if k == "":
            pos_args += f'{v} '
        else:
            ret_str += f'{k} {v} '

    if multi_module:
        ret_str += multiome_suffix

    pos_args+= f'{output[0]} {gene_info_file} '

    # Append positional args at the end of optional args
    pos_args+=ret_str

    return pos_args


def get_mito(wildcards, config):
    if config['chr_prefix'] == None or config['mito_prefix'].startswith(config['chr_prefix']):
        ret_val = config['mito_prefix'] if not config['mito_prefix'].endswith('-') else config['mito_prefix'][:-1]
        return ret_val
    else:
        ret_val = config['chr_prefix'] + config['mito_prefix']
        ret_val =  ret_val if not ret_val.endswith('-') else ret_val[:-1]
        return ret_val


def subset_to_chr(wildcards, config):
    if config['chr_prefix'] == None or str(config['split_bams_pipeline']['subset_chr']).startswith(config['chr_prefix']):
        return config['split_bams_pipeline']['subset_chr']
    else:
        return config['chr_prefix'] + str(config['split_bams_pipeline']['subset_chr'])


def get_mito_file(wildcards, config):
    if config['gt_check']:
        return False
    elif not config['gt_check']:
        return (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['mito_reads_file']}"
            )
    