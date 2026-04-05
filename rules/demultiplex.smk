from lib.io import (
    get_filt_barcodes, 
    get_cellsnp_inputs,
    get_vir_inputs
)
from lib.params import (
    get_params_demux
)


# Need to change this to an elegant solution
def get_inputs_create_h5ad(wildcards):
    ret_list = []
    # Means there's a final count matrix in vireo
    if config['demux_type'].lower() not in [
        'solo', 'vireo', 'both', 
        'add_solo', 'add_vireo'
        ]:
        ret_list.append(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['STARsolo_pipeline']['genefull_matrix']}"
        )

    return ret_list


def get_inputs_demux_solo(wildcards):
    ret_list = []
    # Means there's a final count matrix in vireo
    if config['demux_type'].lower() == 'add_solo':
        ret_list.append(
            f"{config['gt_demux_pipeline']['final_count_matrix_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['gt_demux_pipeline']['final_count_matrix_h5ad']}"
        )
    elif config['demux_type'].lower() == 'solo':
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
    if config['demux_type'].lower() == 'add_vireo':
        ret_list.append(
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"
        )

    elif config['demux_type'].lower() == 'vireo':
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
                    f"{config['fold_struct_gt_demux']}{wildcards.vcf_type}/" # WILDCARDS
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
            f"{config['fold_struct_gt_demux']}{wildcards.vcf_type}/" # WILDCARDS
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



# Resource Allocation ------------------
def allocate_mem_DXP(wildcards, attempt):
    return 3500*attempt+3500


def allocate_time_DXP(wildcards, attempt):
    return 15*attempt+15


# --------------------------------------
rule create_h5ad_only:
    input:
        get_inputs_create_h5ad

    output:
        f"{config['gt_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct']}{config['gt_demux_pipeline']['final_count_matrix_h5ad']}".replace("_vS", ""),
        f"{config['gt_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct']}{config['gt_demux_pipeline']['demultiplex_info']}".replace("_vS", "")

    params:
        pool_name=lambda wildcards: wildcards.pool, # WILDCARDS
        extra=lambda wc, input, output: get_params_demux(input, output, config)

    resources:
        mem_mb=allocate_mem_DXP,
        time_min=allocate_time_DXP

    conda: "../envs/basic_sctools.yaml"

    shell:
        """
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """


rule demux_samples_solo:
    input:
        get_inputs_demux_solo

    output:
        f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}",
        f"{config['hashsolo_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['demultiplex_info']}"

    params:
        pool_name=lambda wildcards: wildcards.pool, # WILDCARDS
        extra=lambda wc, input, output: get_params_demux(input, output, config)

    resources:
        mem_mb=allocate_mem_DXP,
        time_min=allocate_time_DXP

    conda: "../envs/basic_sctools.yaml"
    
    shell: 
        """
        trap 'echo "Received SIGTERM, killing children"; kill 0; wait' SIGTERM
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """


rule demux_samples_vireo:
    input:
        get_inputs_demux_vireo

    output:
        f"{config['gt_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['gt_demux_pipeline']['final_count_matrix_h5ad']}",
        f"{config['gt_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct_demux']}{config['gt_demux_pipeline']['demultiplex_info']}"

    params:
        pool_name=lambda wildcards: wildcards.pool, # WILDCARDS
        extra=lambda wc, input, output: get_params_demux(input, output, config)

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
        pool_name=lambda wildcards: wildcards.pool, # WILDCARDS
        extra=lambda wc, input, output: get_params_demux(input, output, config)

    resources:
        mem_mb=allocate_mem_DXP,
        time_min=allocate_time_DXP

    conda: "../envs/basic_sctools.yaml"
    
    shell: 
        """  
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """
