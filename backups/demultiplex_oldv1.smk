def get_inputs_demux(wildcards):
    ret_list = [
        f"{config['STARsolo_pipeline']['bams_dir']}"
        f"{config['fold_struct']}"
        f"{config['STARsolo_pipeline']['genefull_matrix']}"
        ]
    if global_vars.ONLY_SOLO:
        checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
        ret_list.append(
            f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
            )
    elif global_vars.ONLY_VIREO:
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
        if config['gt_conv'] is not None:
            ret_list.append(config['gt_conv'])

    elif global_vars.BOTH_DEMUX:
        ret_list.append(
            f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
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
        if config['gt_conv'] is not None:
            ret_list.append(config['gt_conv'])
    else:
        raise ValueError(
            "Unexpected inputs to the rule 'demux_samples'! Please "
            "check the INPUTS and the FLAGS!"
            )

    return ret_list


def get_condn(wildcards):
    if global_vars.ONLY_SOLO:
        return 'S'
    elif global_vars.ONLY_VIREO:
        return 'V'
    elif global_vars.BOTH_DEMUX:
        return 'B'
    else:
        return None


def get_inputs_add_demux(wildcards):
    ret_list=[
        f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}"
        f"{config['fold_struct_demux']}"
        f"{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"
        ]
    if global_vars.ADD_SOLO:
        ret_list.append(
            f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}"
            f"{config['fold_struct_demux']}"
            f"{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
            )

    elif global_vars.ADD_VIREO:
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
        if config['gt_conv'] is not None:
            ret_list.append(config['gt_conv'])

    else:
        raise ValueError(
            "Unexpected inputs to the rule "
            "'add_obs_to_final_count_matrix'! Please check the INPUTS "
            "and the FLAGS!"
            )

    return ret_list


def get_condn2(wildcards):
    if global_vars.ADD_SOLO:
        return 'S'
    elif global_vars.ADD_VIREO:
        return 'V'
    else:
        return None


# Resource Allocation ------------------
def allocate_mem_DXP(wildcards, attempt):
    return 3500*attempt+3500


def allocate_time_DXP(wildcards, attempt):
    return 15*attempt+15


def allocate_mem_AOTFCM(wildcards, attempt):
    return 2000+500*(attempt-1)


def allocate_time_AOTFCM(wildcards, attempt):
    return 15+10*(attempt-1)

# --------------------------------------

if global_vars.ONLY_SOLO or global_vars.ONLY_VIREO or global_vars.BOTH_DEMUX:
    rule demux_samples:
        input:
            get_inputs_demux

        output:
            f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}",
            f"{config['hashsolo_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['demultiplex_info']}"

        params:
            mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
            min_genes=config['min_genes_per_cell'], # Min #genes per cell
            min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
            samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
            cols=config['hashsolo_demux_pipeline']['columns_to_pick'],  # Columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
            genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
            pool_name=lambda wildcards: wildcards.id1.replace('-', '_')+'_cDNA',
            hto_sep=config['hashsolo_demux_pipeline']['hto_sep'],
            mito_prefix=config['mito'], # Mitochondrial genes' (names') prefix
            condn=get_condn,
            subid_convert=config['hashsolo_demux_pipeline']['SubID_convert']

        resources:
            mem_mb=allocate_mem_DXP,
            time_min=allocate_time_DXP

        conda: "../envs/basic_sctools.yaml"
        
        shell: 
            """
            read -r -a array <<< "{input}"
            opts_cs=("--calico_solo")
            opts_vs=("--vireo_out" "-converter_file")
            opts_both=(${{opts_cs[@]}} ${{opts_vs[@]}})
            cmd_str=""
            if [ "{params.condn}" == "S" ] && \
               [ "{params.subid_convert}" == "True" ]; then
                if [ "{params.hto_sep}" != "None" ]; then
                    for i in $(seq 1 $((${{#array[@]}}-1)) )
                    do
                        cmd_str+="${{opts_cs[((i-1))]}} ${{array[i]}} "
                    done
                    cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                    python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --demux_info {output[1]} \
                    -m {params.mito} -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
                else
                    python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --demux_info {output[1]} \
                    -m {params.mito} -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
                fi
            elif [ "{params.condn}" == "S" ] && \
                 [ "{params.subid_convert}" != "True" ]; then
                if [ "{params.hto_sep}" != "None" ]; then
                    for i in $(seq 1 $((${{#array[@]}}-1)) )
                    do
                        cmd_str+="${{opts_cs[((i-1))]}} ${{array[i]}} "
                    done
                    cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                    python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --demux_info {output[1]} \
                    -m {params.mito} -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
                else
                    python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --demux_info {output[1]} \
                    -m {params.mito} -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
                fi
            elif [ "{params.condn}" == "V" ]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_vs[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                python3 helper_py_scripts/demul_samples.py \
                {input[0]} {output[0]} {params.genes_info} \
                --demux_info {output[1]} \
                -m {params.mito} -g {params.min_genes} -c {params.min_cells} \
                --mito_prefix {params.mito_prefix} ${{cmd_str}} \
                --pool_name {params.pool_name}

            elif [ "{params.condn}" == "B" ]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_both[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                if [[ "${{#array[@]}}" -lt 4 ]] && [[ "{params.hto_sep}" != "None" ]] && \
                [[ "{params.subid_convert}" == "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
                elif [[ "${{#array[@]}}" -lt 4 ]] && [[ "{params.hto_sep}" != "None" ]] && \
                [[ "{params.subid_convert}" != "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
                elif [[ "${{#array[@]}}" -lt 4 ]] && [[ "{params.hto_sep}" == "None" ]] && \
                [[ "{params.subid_convert}" == "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
                elif [[ "${{#array[@]}}" -lt 4 ]] && [[ "{params.hto_sep}" == "None" ]] && \
                [[ "{params.subid_convert}" != "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
                elif [[ "${{#array[@]}}" -eq 4 ]] && [[ "{params.hto_sep}" != "None" ]] && \
                [[ "{params.subid_convert}" == "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
                elif [[ "${{#array[@]}}" -eq 4 ]] && [[ "{params.hto_sep}" != "None" ]] && \
                [[ "{params.subid_convert}" != "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
                elif [[ "${{#array[@]}}" -eq 4 ]] && [[ "{params.hto_sep}" == "None" ]] && \
                [[ "{params.subid_convert}" == "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
                elif [[ "${{#array[@]}}" -eq 4 ]] && [[ "{params.hto_sep}" == "None" ]] && \
                [[ "{params.subid_convert}" != "True" ]]; then
                    python3 helper_py_scripts/demul_samples.py {input[0]} {output[0]} \
                    {params.genes_info} --demux_info {output[1]} -m {params.mito} \
                    -g {params.min_genes} \
                    -c {params.min_cells} --mito_prefix {params.mito_prefix} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
                fi
            else
                echo "Nothing happened! Check 'condn' in params!"
            fi
            """


# To do the 'ADD_SOLO' part
if global_vars.ADD_VIREO or global_vars.ADD_SOLO:
    # If a previous run of the pipeline has produced final_count_matrix using solo
    # To add results of
    rule add_obs_to_final_count_matrix:
        input:
            get_inputs_add_demux

        output:
            f"{config['gt_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_gt_demux2']}{config['gt_demux_pipeline']['final_count_matrix_h5ad']}"

        resources:
            mem_mb=allocate_mem_AOTFCM,
            time_min=allocate_time_AOTFCM
    
        params:
            samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
            cols=config['hashsolo_demux_pipeline']['columns_to_pick'],  # Columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
            genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
            pool_name=lambda wildcards: wildcards.id1.replace('-', '_')+'_cDNA',
            hto_sep=config['hashsolo_demux_pipeline']['hto_sep'],
            condn=get_condn2,
            subid_convert=config['hashsolo_demux_pipeline']['SubID_convert']

        conda: "../envs/basic_sctools.yaml"

        shell: 
            """
            read -r -a array <<< "{input}"
            opts_cs=("--calico_solo")
            opts_vs=("--vireo_out" "-converter_file")
            cmd_str=""
            if [[ "{params.condn}" == "S" ]] && [[ "{params.hto_sep}" != "None" ]] && \
            [[ "{params.subid_convert}" == "True" ]]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_cs[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --wet_lab_file {params.samples_info} --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
            elif [[ "{params.condn}" == "S" ]] && [[ "{params.hto_sep}" != "None" ]] && \
            [[ "{params.subid_convert}" != "True" ]]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_cs[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --wet_lab_file {params.samples_info} --hto_sep "{params.hto_sep}" \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
            elif [[ "{params.condn}" == "S" ]] && [[ "{params.hto_sep}" == "None" ]] && \
            [[ "{params.subid_convert}" == "True" ]]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_cs[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name}
            elif [[ "{params.condn}" == "S" ]] && [[ "{params.hto_sep}" == "None" ]] && \
            [[ "{params.subid_convert}" != "True" ]]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_cs[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    --wet_lab_file {params.samples_info} \
                    --columns {params.cols} ${{cmd_str}} \
                    --pool_name {params.pool_name} \
                    --no-subid_convert
            elif [[ "{params.condn}" == "V" ]]; then
                for i in $(seq 1 $((${{#array[@]}}-1)) )
                do
                    cmd_str+="${{opts_vs[((i-1))]}} ${{array[i]}} "
                done
                cmd_str="${{cmd_str:0:${{#cmd_str}}-1}}"
                if [[ "${{#array[@]}}" -lt 3 ]]; then
                    python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    ${{cmd_str}} --pool_name {params.pool_name}
                else
                    python3 helper_py_scripts/demul_samples.py \
                    {input[0]} {output[0]} {params.genes_info} \
                    ${{cmd_str}} --pool_name {params.pool_name}
                fi
            else
                echo "Nothing happened! Check 'condn' in params!"
            fi
            """
