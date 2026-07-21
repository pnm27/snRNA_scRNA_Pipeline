#!/usr/bin/env python3

from .utils import read_files_ext
import re, glob2, os

# STARsolo ------------------------------------------------------------------------------
# To check if the last digit of the line in the error log of STARsolo is a number
def check_isnumber(x):
    try:
        int(x)
        return True
   
    except ValueError:
        return False


# We can use this to similarly change other params in the log
def get_limitsjdbval_coll(wildcards, config, resources):
    '''
    This function reads the log file created per attempt to change the parameter "limitsjdbInsertNsj" and "limitOutSJcollapsed" in STARsolo
    Serially produce output as a list in the sequence "limitsjdbInsertNsj", "limitOutSJcollapsed", etc.
    '''
    # This is to check the log file produced after each attempt for the error value
    file_p_temp = f"{config['fold_struct']}".format(**wildcards)
    SS_params_file = "{}{pool}-cDNA.txt".format(config['STARsolo_pipeline']['star_params_dir'], **wildcards)
    
    ins_nsj = 1000000 # DEFAULT
    sj_collap = 1000000 # DEFAULT
    limitbamsortram = 5000000000
    ram_limit_pattern = re.compile(r"--limitBAMsortRAM ([0-9]+)")

    # This is to check the parameters file
    if os.path.isfile(SS_params_file):
        with open(SS_params_file) as fin:
            for line in fin:
                nsj_limit = re.search(r"--limitSjdbInsertNsj (\d+)", line)
                sj_limit = re.search(r"--limitOutSJcollapsed (\d+)", line)
                ram_limit = ram_limit_pattern.search(line)

                temp_nsj = int(nsj_limit.group(1)) if nsj_limit else 0
                temp_sj_coll = int(sj_limit.group(1)) if sj_limit else 0

                ins_nsj = max(
                    temp_nsj,
                    temp_sj_coll,
                    ins_nsj,
                    sj_collap,
                )
                sj_collap = ins_nsj

                if ram_limit:
                    limitbamsortram = max(
                        limitbamsortram,
                        int(ram_limit.group(1)),
                    )


        return {
            'ins_nsj': ins_nsj, 
            'sj_collap': sj_collap, 
            'limitbamsortram': limitbamsortram
        }
    
    log_list = glob2.glob(
        "{}{}_STARsolo_log.txt*".format(config['STARsolo_pipeline']['bams_dir'], 
            file_p_temp)
    )
    # limitbamsortram = int(resources.mem_mb * 0.5 * 1000000) # resources.mem_mb * (resources.cpus_per_task - 1) * 1000000 # DEFAULT
    for log_file in log_list: 
        with open(log_file) as fin:
            for line in fin:
                if line.lower().startswith("solution") and "limitSjdbInsertNsj" in line \
                    and check_isnumber(line.split()[-1]):
                    ins_nsj = max(int(line.split()[-1]), ins_nsj)
                    sj_collap = max(ins_nsj, sj_collap)

                elif line.lower().startswith("solution") and "limitOutSJcollapsed" in line:
                    sj_collap = 1000000*(1+resources.attempt)
                    ins_nsj = max(ins_nsj, sj_collap)
                
                elif line.lower().startswith("solution") and "limitBAMsortRAM" in line:
                    ram_limit = ram_limit_pattern.search(line)
                    if ram_limit:
                        limitbamsortram = max(
                            limitbamsortram,
                            int(ram_limit.group(1)),
                        )
                    
    
    return {
        'ins_nsj': ins_nsj, 
        'sj_collap': sj_collap, 
        'limitbamsortram': limitbamsortram
    }
   

# ---------------------------------------------------------------------------------------

# genotype demux ------------------------------------------------------------------------
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
        

def get_umiTag(wildcards, config):
    # For ATAC return None
    if 'multiome' in config['last_step'].lower() and \
        wildcards.modality.lower() == 'atac':
        return None
    else:
        return config['gt_demux_pipeline']['umi_tag']
    

def get_cmd_str_vireo(wildcards, input, config) -> None | int:
    temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
    # Make snakemake's wildcard same as the value in the "pool" column
    samp_name = wildcards.pool # WILDCARDS
    ret_str = ""

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


# ---------------------------------------------------------------------------------------

# Demultiplex ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------------------

# Split bams ----------------------------------------------------------------------------

def get_params_create_inp_splitBams(wildcards, config):
    ret_str = ''
    donor_params = config['split_bams_pipeline']['donor_name_converter']
    demuxSplit_params = config['split_bams_pipeline']['split_by']

    params_dict = {
        "conv_file": [ donor_params.get('file', None), "--conv_file"],
        "conv_file_FromCol": [ donor_params.get('from_column', None), "--conv_file_from_col"],
        "conv_file_ToCol": [ donor_params.get('to_column', None), "--conv_file_to_col"],
        "demux_split_method": [ demuxSplit_params.get('demux', None), "--split_by"],
        "splitBy_h5adCol": [ demuxSplit_params.get('column', None), "--h5ad_donor_column"],
    }
    
    for _, v in params_dict.items():
        if v[0]:
            ret_str += f'{v[1]} {v[0]} '

    return ret_str


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
    