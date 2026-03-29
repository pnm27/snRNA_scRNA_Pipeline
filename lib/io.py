#!/usr/bin/env python3


import os
from .utils import read_files_ext, ret_cols


# Columns for a vcf_info file when provided with a vcf per pool
col_set = [
    "pool", 
    "n_dons", 
    "vcf"
]
# Make this function to take first set of final_count_matrix file, if one demux is already done
# or to simply filter the count matrix file (output of alignment) according to the same filters
# as in calico_solo's rule
def get_filt_barcodes(wildcards, config, global_vars):

    if global_vars.ADD_SOLO:
        pass
    elif global_vars.ADD_VIREO:
        return (
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"
        )
    elif global_vars.ONLY_VIREO or global_vars.BOTH_DEMUX:
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
    ret_list=[
        f"{config['gt_demux_pipeline']['inp_for_cellsnp_dir']}"
        f"{config['fold_struct_filt_bc']}.txt"
    ]
    if 'multiome' in config['last_step'].lower():
        if 'cdna' in wildcards.pool.lower(): # WILDCARDS
            ret_list.append((
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
                ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            )
        else:
            ret_list.append((
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
                ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            )
    else:
        ret_list.append((
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"
        ))
    
    # 'cellsnp_ref_snps' is useful when per-Pool vcfs need to subset for SNVs
    if config['gt_demux_pipeline']['cellsnp_ref_snps'] is not None:
        # modify this to accept either a text file containing per pool vcf or one single vcf file
        if config['gt_demux_pipeline']['cellsnp_ref_snps'].endswith('.vcf.gz'):
            ret_list.append(config['gt_demux_pipeline']['cellsnp_ref_snps'])
            ret_list.append(config['gt_demux_pipeline']['genome_1k_ref'])
        
        # TEMPORARY: When vireo has per-pool runs and we want similar pileups
        else:
            col_set1 = ["pool", "n_dons", "vcf"]
            n_cols = ret_cols(config['gt_demux_pipeline']['vcf_info'])
            samp_name = '-'.join(wildcards.pool.split('-')[:-1])
            if config['gt_demux_pipeline']['donors_per_pool_file'] is None \
            and n_cols == 3:
                temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'],
                          names=col_set1)
                ret_list.append(
                    temp_df.loc[temp_df["pool"].str.contains(samp_name,na=False), "vcf"]
                    .values[0]
                )

    elif config['gt_demux_pipeline']['cellsnp_ref_snps'] is None:
        # If no ref provided then use 1000 Genomes Project vcf as input
        ret_list.append(config['gt_demux_pipeline']['genome_1k_ref'])


    return ret_list

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
    ret_list = [
        f"{config['gt_demux_pipeline']['cellsnp_dir']}{config['fold_struct_gt_demux']}"
        f"{config['gt_demux_pipeline']['cellsnp_cells']}"
        ]
    # if config['gt_demux_pipeline']['vcf_info'].endswith('vcf.gz'):
    #     ret_dict['donor_info'] = config['gt_demux_pipeline']['vcf_info']
    # elif config['gt_demux_pipeline']['vcf_info'].endswith('vcf'):
    #     raise IOError("Given vcf file in 'vcf_info' in config file is not BGZF formatted!")
    # elif config['gt_demux_pipeline']['vcf_info'].endswith('txt') or config['gt_demux_pipeline']['vcf_info'].endswith('csv') \
    #     or config['gt_demux_pipeline']['vcf_info'].endswith('tsv'):
    if config['gt_demux_pipeline']['vcf_info'] is not None \
        and os.path.isfile(config['gt_demux_pipeline']['vcf_info']):
        n_cols = ret_cols(config['gt_demux_pipeline']['vcf_info'])
        # Make snakemake's wildcard same as the value in the "pool" column
        samp_name = '-'.join(wildcards.pool.split('-')[:-1]) # WILDCARDS
        # For condn. 1)
        if config['gt_demux_pipeline']['donors_per_pool_file'] is None \
            and n_cols == 3:
            temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'], 
                                     names=col_set)
            try:
                ret_list.append(
                    temp_df.loc[temp_df["pool"] == samp_name, "vcf"]
                    .values[0]
                    )
            except:
                pass
        # Run without vcf
        elif config['gt_demux_pipeline']['donors_per_pool_file'] is None \
            and n_cols == 2:
            pass
        # For condn. 2)
        elif config['gt_demux_pipeline']['donors_per_pool_file'].endswith('vcf.gz') \
            and n_cols == 3:
            raise ValueError("Expecting 2-columned"\
                " file: pool name and expected donors")
        elif config['gt_demux_pipeline']['donors_per_pool_file'].endswith('vcf.gz') \
            and n_cols == 2:
            ret_list.append(config['gt_demux_pipeline']['donors_per_pool_file'])
        elif config['gt_demux_pipeline']['donors_per_pool_file'].endswith('vcf'):
            raise IOError("Given vcf file in 'donors_per_pool_file' in "
                          "config file is not BGZF formatted!")
        else:
            raise IOError("Unexpected file type or value provided to "
                          "'donors_per_pool_file' in config file")

    return ret_list

