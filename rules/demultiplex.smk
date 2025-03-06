def get_inputs_demux_solo(wildcards):
    ret_list = []
    # Means there's a final count matrix in vireo
    if global_vars.ADD_SOLO:
        ret_list.append(
            f"{config['gt_demux_pipeline']['final_count_matrix_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['gt_demux_pipeline']['final_count_matrix_h5ad']}"
        )
    elif global_vars.ONLY_SOLO:
        ret_list.append(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['STARsolo_pipeline']['genefull_matrix']}"
        )

    ret_list.append(
        f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}"
        f"{config['fold_struct_demux']}"
        f"{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
        )
    
    return ret_list


def get_inputs_demux_vireo(wildcards):
    multi_module = ['multiome', 'gt_demux']
    ret_list = []
    # Means there's a final count matrix in solo
    if global_vars.ADD_VIREO:
        ret_list.append(
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"
        )

    elif global_vars.ONLY_VIREO:
        if all([ m in config['last_step'].lower() for m in multi_module]):
            ret_list.extend([
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/filtered_feature_bc_matrix/matrix.mtx.gz",
                f"{config['gt_demux_pipeline']['vireosnp_dir']}"
                f"{config['fold_struct_gt_demux']}ATAC/"
                f"{config['gt_demux_pipeline']['donors_classification']}",
                f"{config['gt_demux_pipeline']['vireosnp_dir']}"
                f"{config['fold_struct_gt_demux']}cDNA/"
                f"{config['gt_demux_pipeline']['donors_classification']}"
            ])
        else:
            ret_list.append(
                f"{config['STARsolo_pipeline']['bams_dir']}"
                f"{config['fold_struct']}"
                f"{config['STARsolo_pipeline']['genefull_matrix']}"
            )
            if config['last_step'].lower().endswith('multi_vcf'):
                ret_list.append(
                    f"{config['gt_demux_pipeline']['vireosnp_dir']}"
                    f"{config['fold_struct_gt_demux']}{wildcards.vcf_type}/"
                    f"{config['gt_demux_pipeline']['donors_classification']}"
                    )
            else:
                ret_list.append(
                    f"{config['gt_demux_pipeline']['vireosnp_dir']}"
                    f"{config['fold_struct_gt_demux']}"
                    f"{config['gt_demux_pipeline']['donors_classification']}"
                    )
    if config['gt_demux_pipeline']['donorName_conv']['file'] is not None:
        ret_list.append(config['gt_demux_pipeline']['donorName_conv']['file'])

    return ret_list

def get_inputs_demux_both(wildcards):
    
    ret_list = [
        f"{config['STARsolo_pipeline']['bams_dir']}"
        f"{config['fold_struct']}"
        f"{config['STARsolo_pipeline']['genefull_matrix']}",
        f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}"
        f"{config['fold_struct_demux']}"
        f"{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
        ]

    if config['last_step'].lower().endswith('multi_vcf')  :
        ret_list.append(
            f"{config['gt_demux_pipeline']['vireosnp_dir']}"
            f"{config['fold_struct_gt_demux']}{wildcards.vcf_type}/"
            f"{config['gt_demux_pipeline']['donors_classification']}"
            )
    else:
        ret_list.append(
            f"{config['gt_demux_pipeline']['vireosnp_dir']}"
            f"{config['fold_struct_gt_demux']}"
            f"{config['gt_demux_pipeline']['donors_classification']}"
            )
    if config['gt_demux_pipeline']['donorName_conv']['file'] is not None:
        ret_list.append(config['gt_demux_pipeline']['donorName_conv']['file'])


    return ret_list


def get_params(wildcards, input, output):
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
    }
    gene_info_file = config['gene_info_file']
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

    if global_vars.ONLY_SOLO or global_vars.ADD_SOLO:
        inp_files += solo_inp
    elif global_vars.ONLY_VIREO or global_vars.ADD_VIREO:
        if multi_module:
            inp_files += multiome_inp
        else:
            inp_files += vireo_inp
    elif global_vars.BOTH_DEMUX:
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


# DEPRACATED
# def get_inputs_add_demux(wildcards):
#     ret_list=[
#         f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}"
#         f"{config['fold_struct_demux']}"
#         f"{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"
#         ]
#     if global_vars.ADD_SOLO:
#         ret_list.append(
#             f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}"
#             f"{config['fold_struct_demux']}"
#             f"{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
#             )

#     elif global_vars.ADD_VIREO:
#         if config['last_step'].lower().endswith('multi_vcf'):
#             ret_list.append(
#                 f"{config['gt_demux_pipeline']['vireosnp_dir']}"
#                 f"{config['fold_struct_gt_demux']}{wildcards.vcf_type}/"
#                 f"{config['gt_demux_pipeline']['donors_classification']}"
#                 )
#         else:
#             ret_list.append(
#                 f"{config['gt_demux_pipeline']['vireosnp_dir']}"
#                 f"{config['fold_struct_gt_demux']}"
#                 f"{config['gt_demux_pipeline']['donors_classification']}"
#                 )
#         if config['gt_demux_pipeline']['donorName_conv']['file'] is not None:
#             ret_list.append(config['gt_demux_pipeline']['donorName_conv']['file'])

#     else:
#         raise ValueError(
#             "Unexpected inputs to the rule "
#             "'add_obs_to_final_count_matrix'! Please check the INPUTS "
#             "and the FLAGS!"
#             )

#     return ret_list


# DEPRACATED
# def get_condn2(wildcards):
#     if global_vars.ADD_SOLO:
#         return 'S'
#     elif global_vars.ADD_VIREO:
#         return 'V'
#     else:
#         return None


# Resource Allocation ------------------
def allocate_mem_DXP(wildcards, attempt):
    return 3500*attempt+3500


def allocate_time_DXP(wildcards, attempt):
    return 15*attempt+15

# DEPRACATED
# def allocate_mem_AOTFCM(wildcards, attempt):
#     return 2000+500*(attempt-1)

# DEPRACATED
# def allocate_time_AOTFCM(wildcards, attempt):
#     return 15+10*(attempt-1)

# --------------------------------------

# if global_vars.ONLY_SOLO or global_vars.ONLY_VIREO or global_vars.BOTH_DEMUX:
rule demux_samples_solo:
    input:
        get_inputs_demux_solo

    output:
        f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}",
        f"{config['hashsolo_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['demultiplex_info']}"

    params:
        # mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        # min_genes=config['min_genes_per_cell'], # Min #genes per cell
        # min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        # samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
        # cols=config['hashsolo_demux_pipeline']['columns_to_pick'],  # Columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
        # genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
        pool_name=lambda wildcards: wildcards.pool.replace('-', '_')+'_cDNA', #WILDCARDS
        # hto_sep=config['hashsolo_demux_pipeline']['hto_sep'],
        # mito_prefix=config['mito_prefix'], # Mitochondrial genes' (names') prefix
        # condn=get_condn,
        # subid_convert=config['hashsolo_demux_pipeline']['SubID_convert']
        extra=get_params

    resources:
        mem_mb=allocate_mem_DXP,
        time_min=allocate_time_DXP

    conda: "../envs/basic_sctools.yaml"
    
    shell: 
        """  
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """


rule demux_samples_vireo:
    input:
        get_inputs_demux_vireo

    output:
        f"{config['gt_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['gt_demux_pipeline']['final_count_matrix_h5ad']}",
        f"{config['gt_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct_demux']}{config['gt_demux_pipeline']['demultiplex_info']}"

    params:
        # mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        # min_genes=config['min_genes_per_cell'], # Min #genes per cell
        # min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        # samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
        # cols=config['hashsolo_demux_pipeline']['columns_to_pick'],  # Columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
        # genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
        pool_name=lambda wildcards: wildcards.pool.replace('-', '_')+'_cDNA', #WILDCARDS
        # hto_sep=config['hashsolo_demux_pipeline']['hto_sep'],
        # mito_prefix=config['mito_prefix'], # Mitochondrial genes' (names') prefix
        # condn=get_condn,
        # subid_convert=config['hashsolo_demux_pipeline']['SubID_convert']
        extra=get_params

    resources:
        mem_mb=allocate_mem_DXP,
        time_min=allocate_time_DXP

    conda: "../envs/basic_sctools.yaml"
    
    shell: 
        """  
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """


rule demux_samples_both:
    input:
        get_inputs_demux_both

    output:
        f"{config['demultiplex']['demux_count_matrix_dir']}{config['fold_struct_demux']}{config['demultiplex']['final_count_matrix_h5ad']}",
        f"{config['demultiplex']['demux_count_matrix_dir']}{config['fold_struct_demux']}{config['demultiplex']['demultiplex_info']}"

    params:
        # mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        # min_genes=config['min_genes_per_cell'], # Min #genes per cell
        # min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        # samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
        # cols=config['hashsolo_demux_pipeline']['columns_to_pick'],  # Columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
        # genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
        pool_name=lambda wildcards: wildcards.pool.replace('-', '_')+'_cDNA', #WILDCARDS
        # hto_sep=config['hashsolo_demux_pipeline']['hto_sep'],
        # mito_prefix=config['mito_prefix'], # Mitochondrial genes' (names') prefix
        # condn=get_condn,
        # subid_convert=config['hashsolo_demux_pipeline']['SubID_convert']
        extra=get_params

    resources:
        mem_mb=allocate_mem_DXP,
        time_min=allocate_time_DXP

    conda: "../envs/basic_sctools.yaml"
    
    shell: 
        """  
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """



# DEPRACATED
# if global_vars.ADD_VIREO or global_vars.ADD_SOLO:
#     # If a previous run of the pipeline has produced final_count_matrix using solo
#     # To add results of
#     rule add_obs_to_final_count_matrix:
#         input:
#             get_inputs_add_demux

#         output:
#             f"{config['gt_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_gt_demux2']}{config['gt_demux_pipeline']['final_count_matrix_h5ad']}"

#         resources:
#             mem_mb=allocate_mem_AOTFCM,
#             time_min=allocate_time_AOTFCM
    
#         params:
#             # samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
#             # cols=config['hashsolo_demux_pipeline']['columns_to_pick'],  # Columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
#             # genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
#             pool_name=lambda wildcards: wildcards.pool.replace('-', '_')+'_cDNA', #WILDCARDS
#             # hto_sep=config['hashsolo_demux_pipeline']['hto_sep'],
#             # condn=get_condn2,
#             # subid_convert=config['hashsolo_demux_pipeline']['SubID_convert']
#             extra=get_params

#         conda: "../envs/basic_sctools.yaml"

#         shell: 
#             """
#             python3 helper_py_scripts/demul_samples.py {params.extra}
#             """
