import glob2, re
import pandas as pd


def get_r1_hto_fastqs(wildcards):
    temp_fold=f"{config['fold_struct']}".format(id1=wildcards.id1).replace('-cDNA', '-HTO')
    r1_files = sorted((
        glob2.glob("{parent_dir}{fs}*{suff}"
        .format(parent_dir=config['HTO_fastqs_dir'], fs=temp_fold, suff=config['R1_suffix']))
        ))
    return r1_files


def get_r2_hto_fastqs(wildcards):
    temp_fold=f"{config['fold_struct']}".format(id1=wildcards.id1).replace('-cDNA', '-HTO')
    r2_files = sorted((
        glob2.glob("{parent_dir}{fs}*{suff}"
        .format(parent_dir=config['HTO_fastqs_dir'], fs=temp_fold, suff=config['R2_suffix']))
        ))
    return r2_files



def get_hto_fastqs(wildcards):
    temp_fold=f"{config['fold_struct']}".format(id1=wildcards.id1).replace('-cDNA', '-HTO')
    all_files = sorted((
        glob2.glob("{parent_dir}{fs}*{suff}"
        .format(parent_dir=config['HTO_fastqs_dir'], fs=temp_fold, suff=config['R1_suffix'].replace('1', '*')))
        ))
    return all_files


rule create_FB:
    input:
        R1=get_r1_hto_fastqs,
        R2=get_r2_hto_fastqs,
        sample_sheet=config['wet_lab_info']
        
    # priority: 10

    output:
        f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['feature_barcodes']}"

    params:
        sample_name=lambda wildcards: wildcards.id1.replace('-', '_') + '_cDNA'

    resources:
        mem_mb=allocate_mem_CFB, #allocate_mem_KBP,
        time_min=allocate_time_CFB

    shell: 
        """
        python3 helper_py_scripts/create_Feat_Barc.py {input.sample_sheet} \
        -o {output} -s {params.sample_name} -c {config[kb_pipeline][columns_to_pick]}
        """


         

if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':
    rule create_mismatch_fasta:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['feature_barcodes']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['features_mismatch_fa']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['features_mismatch_t2g']}"

        params:
            headers=config['kb_pipeline']['headers'] # Does the feature barcodes file hasve headers

        resources:
            mem_mb=allocate_mem_CMF, #allocate_mem_KBP,
            time_min=allocate_time_CMF


        shell: 
            """
            if [[ "{params.headers}" == "yes" ]] || [[ "{params.headers}" == "True" ]]; then
                python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]} --header

            else
                python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]}

            fi
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':
    rule create_mismatch_fasta:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['feature_barcodes']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['features_mismatch_fa']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['features_mismatch_t2g']}"

        params:
            headers=config['kb_pipeline']['headers'] # Does the feature barcodes file hasve headers

        resources:
            mem_mb=allocate_mem_CMF, #allocate_mem_KBP,
            time_min=allocate_time_CMF


        shell: 
            """
            if [[ "{params.headers}" == "yes" ]] || [[ "{params.headers}" == "True" ]]; then
                python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]} --header

            else
                python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]}

            fi
            """

if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':
    rule build_kallisto_index:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['features_mismatch_fa']}"

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['features_mismatch_idx']}"

        params:
            k_mer_len=config['kb_pipeline']['k_mer_length']

        resources:
            mem_mb=allocate_mem_BKI, #allocate_mem_KBP,
            time_min=allocate_time_BKI

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            ml kallisto
            kallisto index -i {output} -k {params.k_mer_len} {input}
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':
    rule build_kallisto_index:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['features_mismatch_fa']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['features_mismatch_idx']}"

        params:
            k_mer_len=config['kb_pipeline']['k_mer_length']

        resources:
            mem_mb=allocate_mem_BKI, #allocate_mem_KBP,
            time_min=allocate_time_BKI

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            ml kallisto
            kallisto index -i {output} -k {params.k_mer_len} {input}
            """


if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':
    rule run_kallisto:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['features_mismatch_idx']}",
            fastq_files=get_hto_fastqs

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_file']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['ec_matrix']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['tx']}"

        params:
            output_pref= lambda wildcards, output: output[0].replace(f"{config['kb_pipeline']['bus_file']}", ''),
            chemistry=config['kb_pipeline']['chemistry']

        threads: 1

        resources:
            mem_mb=allocate_mem_RK, #allocate_mem_KBP,
            time_min=allocate_time_RK

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            ml kallisto
            kallisto bus -i {input[0]} -o {params.output_pref} -x {params.chemistry} -t 8 {input.fastq_files}
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':
    rule run_kallisto:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['features_mismatch_idx']}",
            fastq_files=get_hto_fastqs

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_file']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['ec_matrix']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['tx']}"

        params:
            output_pref= lambda wildcards, output: output[0].replace(f"{config['kb_pipeline']['bus_file']}", ''),
            chemistry=config['kb_pipeline']['chemistry']

        threads: 1

        resources:
            mem_mb=allocate_mem_RK, #allocate_mem_KBP,
            time_min=allocate_time_RK

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            ml kallisto
            kallisto bus -i {input[0]} -o {params.output_pref} -x {params.chemistry} -t 8 {input.fastq_files}
            """


if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':
    rule run_bustools_correct:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_file']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_file_corrected']}"

        params:
            whitelist=config['whitelist']

        resources:
            mem_mb=allocate_mem_RBCor, #allocate_mem_KBP,
            time_min=allocate_time_RBCor

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            ml bustools
            bustools correct -w {params.whitelist} {input} -o {output}
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':
    rule run_bustools_correct:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_file']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_file_corrected']}"

        params:
            whitelist=config['whitelist']

        resources:
            mem_mb=allocate_mem_RBCor, #allocate_mem_KBP,
            time_min=allocate_time_RBCor

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            ml bustools
            bustools correct -w {params.whitelist} {input} -o {output}
            """


if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':        
    rule run_bustools_sort:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_file_corrected']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_file_sorted']}"

        threads: 1

        resources:
            mem_mb=allocate_mem_RBS, #allocate_mem_KBP,
            time_min=allocate_time_RBS

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            ml bustools
            bustools sort -t 4 -o {output} {input}
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':        
    rule run_bustools_sort:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_file_corrected']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_file_sorted']}"

        threads: 1

        resources:
            mem_mb=allocate_mem_RBS, #allocate_mem_KBP,
            time_min=allocate_time_RBS

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            ml bustools
            bustools sort -t 4 -o {output} {input}
            """


if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':
    rule run_bustools_count:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_file_sorted']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['ec_matrix']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['tx']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['features_mismatch_t2g']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_mtx']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_genes']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_barcodes']}"

        params:
            output_pref=lambda wildcards, output: output[0].replace(f"{config['kb_pipeline']['bus_count_mtx']}", '')

        resources:
            mem_mb=allocate_mem_RBCnt, #allocate_mem_KBP,
            time_min=allocate_time_RBCnt

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            ml bustools
            bustools count -o {params.output_pref} --genecounts -g {input[3]} -e {input[1]} -t {input[2]} {input[0]}
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':
    rule run_bustools_count:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_file_sorted']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['ec_matrix']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['tx']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['features_mismatch_t2g']}"

        # priority: 10

        output:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_mtx']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_genes']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_barcodes']}"

        params:
            output_pref=lambda wildcards, output: output[0].replace(f"{config['kb_pipeline']['bus_count_mtx']}", '')

        resources:
            mem_mb=allocate_mem_RBCnt, #allocate_mem_KBP,
            time_min=allocate_time_RBCnt

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            ml bustools
            bustools count -o {params.output_pref} --genecounts -g {input[3]} -e {input[1]} -t {input[2]} {input[0]}
            """


if config['hto_demux_type'] is not None and config['hto_demux_type'].lower() == 'single':
    rule create_h5ad_bustools:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_mtx']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_genes']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_barcodes']}"

        # priority: 10

        output:
            f"{config['demux_pipeline']['h5ad_bustools_dir']}{config['fold_struct_demux']}{config['demux_pipeline']['bustools_h5ad']}"

        resources:
            mem_mb=allocate_mem_CHB, #allocate_mem_KBP,
            time_min=allocate_time_CHB

        shell:
            """
            python3 helper_py_scripts/create_h5ad_from_bustools.py {input[0]} {input[1]} {input[2]} -o {output}
            """
elif config['hto_demux_type'] is not None and config['hto_demux_type'].lower() != 'single':
    rule create_h5ad_bustools:
        input:
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_mtx']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_genes']}",
            f"{config['kb_pipeline']['kallisto_bustools_dir']}{config['fold_struct_kb']}{{i}}{config['kb_pipeline']['bus_count_dir']}{config['kb_pipeline']['bus_count_barcodes']}"

        # priority: 10

        output:
            f"{config['demux_pipeline']['h5ad_bustools_dir']}{config['fold_struct_demux']}_{{i}}{config['demux_pipeline']['bustools_h5ad']}"

        resources:
            mem_mb=allocate_mem_CHB, #allocate_mem_KBP,
            time_min=allocate_time_CHB

        conda: "../envs/basic_sctools.yaml"
        
        shell:
            """
            python3 helper_py_scripts/create_h5ad_from_bustools.py {input[0]} {input[1]} {input[2]} -o {output}
            """
