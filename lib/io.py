#!/usr/bin/env python3

import os
from .utils import read_files_ext, ret_cols


# Columns for a vcf_info file when provided with a vcf per pool
col_set = [
    "pool", 
    "n_dons", 
    "donors", # NEW ADDITION (comma-sep donor list)
    "vcf",
]
# Make this function to take first set of final_count_matrix file, if one demux is already done
# or to simply filter the count matrix file (output of alignment) according to the same filters
# as in calico_solo's rule
def get_filt_barcodes(wildcards, config):

    if config['demux_type'].lower() == 'add_solo':
        pass
    elif config['demux_type'].lower() == 'add_vireo':
        return (
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"
        )
    elif config['demux_type'].lower() == 'vireo' or config['demux_type'].lower() == 'both':
        if 'multiome' in config['last_step'].lower():
            return (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}"
                "/filtered_feature_bc_matrix/matrix.mtx.gz"
            ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS

        else:
            return (
                f"{config['STARsolo_pipeline']['bams_dir']}"
                f"{config['fold_struct']}"
                f"{config['STARsolo_pipeline']['genefull_matrix']}"
            )
    else:
        raise ValueError("Unexpected inputs to the rule 'demux_samples_calico_solo_STARsolo'! Please check the INPUTS and the FLAGS!")
    

def get_cellsnp_inputs(wildcards, config):
    ret_dict={}
    ret_dict['barcodesFile'] = (
        f"{config['gt_demux_pipeline']['inp_for_cellsnp_dir']}"
        f"{config['fold_struct_filt_bc']}.txt"
    )
    if 'multiome' in config['last_step'].lower():
        if 'cdna' in wildcards.pool.lower(): # WILDCARDS
            ret_dict['bam'] = ((
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
                ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            )
        else:
            ret_dict['bam'] = ((
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
                ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            )
    else:
        ret_dict.update({
            'bam': (
                f"{config['STARsolo_pipeline']['bams_dir']}"
                f"{config['fold_struct']}"
                f"{config['STARsolo_pipeline']['bam']}"
            ),
            'bai': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['STARsolo_pipeline']['bai']}"
            )
        })
    
    # 'cellsnp_ref_snps' is useful when per-Pool vcfs need to subset for SNVs
    if config['gt_demux_pipeline']['vcf_info'] is not None \
        and os.path.isfile(config['gt_demux_pipeline']['vcf_info']):
            samp_name = '-'.join(wildcards.pool)
            temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
            n_cols = temp_df.shape[1]
            # For demux run wo_gt, we generally won't have per-donor vcf
            if config['demux_run_type'] == 'wo_gt':
                ret_dict['regionsFile'] = config['gt_demux_pipeline']['genome_1k_ref']
            # If there are more than 3 columns then columns 4 and above
            # have per-donor vcfs 
            # (the 'above' is useful for multi-vcf setup) - TODO
            elif n_cols == 4:
                temp_df = temp_df.iloc[:, :4]
                temp_df.columns = col_set
                ret_dict['regionsFile'] = temp_df.loc[
                    temp_df["pool"] == samp_name, "vcf"].values[0]

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
    ret_dict['cellsnpCells'] = (
        f"{config['gt_demux_pipeline']['cellsnp_dir']}"
        f"{config['fold_struct_gt_demux']}"
        f"{config['gt_demux_pipeline']['cellsnp_cells']}"    
    )

    if config['gt_demux_pipeline']['vcf_info'] is not None \
        and os.path.isfile(config['gt_demux_pipeline']['vcf_info']):
        temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
        n_cols = temp_df.shape[1]
        # Make snakemake's wildcard same as the value in the "pool" column
        samp_name = wildcards.pool # WILDCARDS
        # For condn. 1)
        if n_cols == 4:
            temp_df = temp_df.iloc[:, :n_cols]
            temp_df.columns = col_set
            ret_dict['donorFile'] = (
                temp_df.loc[temp_df["pool"] == samp_name, "vcf"]
                .values[0]
                )
        else:
            ret_dict['donorFile'] = []  

    return ret_dict



def get_inp_splitBam(wildcards, config):
    if config['split_bams_pipeline']['split_by']['input'].lower() == 'raw':
        if config['split_bams_pipeline']['split_by']['demux'].lower() in ['vireo', 'vs']:
            return f"{config['gt_demux_pipeline']['vireosnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['donors_classification']}"
        elif config['split_bams_pipeline']['split_by']['demux'].lower() in \
            ["cs", "calico", "calico_solo", "hashsolo"]:
            return f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
    else:
        # if len(config['split_bams_pipeline']['split_by']['demux']) == 1:
        #     if config['split_bams_pipeline']['split_by']['demux'][0].lower() in ['vireo', 'vs'] or \
        #       config['split_bams_pipeline']['split_by']['demux'][0].lower() in \
        #       ["cs", "calico", "calico_solo", "hashsolo"]:
        #         return f"{config['demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['demux_pipeline']['final_count_matrix_h5ad']}"
        # else:
        return f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"


def get_bam(wildcards, config):
    if 'multiome' in config['last_step'].lower():
        if 'cdna' in wildcards.pool.lower(): # WILDCARDS
            return {
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
                ).format(pool=wildcards.pool.split('/')[0]), # WILDCARDS
                'bai': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bai']}"
                ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            } 
                
        else:
            return {
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
                ).format(pool=wildcards.pool.split('/')[0]), # WILDCARDS
                'bai': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bai']}"
                ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            } 
    else:
        return {
            'bam': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['STARsolo_pipeline']['bam']}"
            ),
            'bai': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['STARsolo_pipeline']['bai']}"
            )
        }
    

def get_bam_to_split(wildcards, config):
    # fullBam = ""
    # subBam = ""
    # filtBam = ""
    fullBam = {}
    subBam = {}
    filtBam = {}
    if 'multiome' in config['last_step'].lower():
        if 'cdna' in wildcards.pool.lower(): # WILDCARDS
            fullBam.update({
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
                ).format(pool=wildcards.pool.split('/')[0])}
            )
            subBam.update({
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
                ).format(
                    pool=wildcards.pool.split('/')[0]).replace(  # WILDCARDS
                        '.bam', config['split_bams_pipeline']['short_bam'])}
            ) 
            filtBam.update({
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
                ).format(
                    pool=wildcards.pool.split('/')[0]).replace( # WILDCARDS
                        '.bam', config['split_bams_pipeline']['filt_bam'])}
            )
        else:
            fullBam.update({
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
                ).format(pool=wildcards.pool.split('/')[0])} # WILDCARDS
            ) 
            subBam.update({
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
                ).format(
                    pool=wildcards.pool.split('/')[0]).replace( # WILDCARDS
                        '.bam', config['split_bams_pipeline']['short_bam'])}
            )
            filtBam.update({
                'bam': (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
                ).format(
                    pool=wildcards.pool.split('/')[0]).replace( # WILDCARDS
                        '.bam', config['split_bams_pipeline']['filt_bam'])}
            )
    else:
        fullBam.update({
            'bam':(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"
            ),
            'bai':(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['STARsolo_pipeline']['bai']}"
            )
            }
        )
        subBam.update({
            'bam': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['split_bams_pipeline']['short_bam']}"
            ),
            'bai': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['split_bams_pipeline']['short_bam']}.bai"
            )
            }
        )
        filtBam.update({
            'bam': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}"
            ),
            'bai': (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}.bai"
            )
            }
        )
        
    if config['gt_check']:
        if config['split_bams_pipeline']['subset_chr'] is None:
            return fullBam
        else:
            return subBam
    else:
        return filtBam