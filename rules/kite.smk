import glob2, re
import pandas as pd


# HTO multiplexing: one HTO per cell ('single') or more ('multiple')
from lib.utils import run_settings
SINGLE_HTO = run_settings(config).hto_multiplexing == 'single'
MULTI_HTO = run_settings(config).hto_multiplexing == 'multiple'


def get_r1_hto_fastqs(wildcards):
    temp_fold=f"{config['modules']['STARsolo']['fold_struct']}".format(**wildcards).replace('-cDNA', '-HTO')
    r1_files = sorted((
        glob2.glob("{parent_dir}{fs}*{suff}"
        .format(parent_dir=config['experiment']['inputs']['fastqs_dir']['HTO'], fs=temp_fold, suff=config['experiment']['inputs']['R1_suffix']))
        ))
    return r1_files


def get_r2_hto_fastqs(wildcards):
    temp_fold=f"{config['modules']['STARsolo']['fold_struct']}".format(**wildcards).replace('-cDNA', '-HTO')
    r2_files = sorted((
        glob2.glob("{parent_dir}{fs}*{suff}"
        .format(parent_dir=config['experiment']['inputs']['fastqs_dir']['HTO'], fs=temp_fold, suff=config['experiment']['inputs']['R2_suffix']))
        ))
    return r2_files



def get_hto_fastqs(wildcards):
    temp_fold=f"{config['modules']['STARsolo']['fold_struct']}".format(**wildcards).replace('-cDNA', '-HTO')
    all_files = sorted((
        glob2.glob("{parent_dir}{fs}*{suff}"
        .format(parent_dir=config['experiment']['inputs']['fastqs_dir']['HTO'], fs=temp_fold, suff=config['experiment']['inputs']['R1_suffix'].replace('1', '*')))
        ))
    return all_files


rule create_FB:
    input:
        R1=get_r1_hto_fastqs,
        R2=get_r2_hto_fastqs
        
    # priority: 10

    output:
        f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['feature_barcodes']}"

    params:
        sample_name=lambda wildcards: wildcards.pool.replace('-', '_') + '_cDNA', # WILDCARDS
        sample_sheet=config['experiment']['inputs']['wet_lab_info']

    resources:
        **get_resources("create_FB")

    shell: 
        """
        python3 helper_py_scripts/create_Feat_Barc.py {params.sample_sheet} \
        -o {output} -s {params.sample_name} -c {config[modules][kite][columns_to_pick]}
        """
  

if SINGLE_HTO:
    rule create_mismatch_fasta:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['feature_barcodes']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['features_mismatch_fa']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['features_mismatch_t2g']}"

        params:
            headers=config['modules']['kite']['headers'] # Does the feature barcodes file hasve headers

        # For snakemake < v8
        threads: 1

        resources:
            **get_resources("create_mismatch_fasta")


        shell: 
            """
            if [[ "{params.headers}" == "yes" ]] || [[ "{params.headers}" == "True" ]]; then
                python3 {config[references][tools][featuremap_script]} {input} --t2g {output[1]} --fa {output[0]} --header

            else
                python3 {config[references][tools][featuremap_script]} {input} --t2g {output[1]} --fa {output[0]}

            fi
            """

elif MULTI_HTO:
    rule create_mismatch_fasta:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['feature_barcodes']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['features_mismatch_fa']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['features_mismatch_t2g']}"

        params:
            headers=config['modules']['kite']['headers'] # Does the feature barcodes file hasve headers

        resources:
            **get_resources("create_mismatch_fasta")


        shell: 
            """
            if [[ "{params.headers}" == "yes" ]] || [[ "{params.headers}" == "True" ]]; then
                python3 {config[references][tools][featuremap_script]} {input} --t2g {output[1]} --fa {output[0]} --header

            else
                python3 {config[references][tools][featuremap_script]} {input} --t2g {output[1]} --fa {output[0]}

            fi
            """

if SINGLE_HTO:
    rule build_kallisto_index:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['features_mismatch_fa']}"

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['features_mismatch_idx']}"

        params:
            k_mer_len=config['modules']['kite']['k_mer_length']

        resources:
            **get_resources("build_kallisto_index")

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            kallisto index -i {output} -k {params.k_mer_len} {input}
            """

elif MULTI_HTO:
    rule build_kallisto_index:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['features_mismatch_fa']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['features_mismatch_idx']}"

        params:
            k_mer_len=config['modules']['kite']['k_mer_length']

        resources:
            **get_resources("build_kallisto_index")

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            kallisto index -i {output} -k {params.k_mer_len} {input}
            """


if SINGLE_HTO:
    rule run_kallisto:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['features_mismatch_idx']}",
            fastq_files=get_hto_fastqs

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_file']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['ec_matrix']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['tx']}"

        params:
            output_pref= lambda wildcards, output: output[0].replace(f"{config['modules']['kite']['bus_file']}", ''),
            chemistry=config['modules']['kite']['chemistry']

        resources:
            **get_resources("run_kallisto")

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            kallisto bus -i {input[0]} -o {params.output_pref} -x {params.chemistry} -t 8 {input.fastq_files}
            """

elif MULTI_HTO:
    rule run_kallisto:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['features_mismatch_idx']}",
            fastq_files=get_hto_fastqs

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_file']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['ec_matrix']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['tx']}"

        params:
            output_pref= lambda wildcards, output: output[0].replace(f"{config['modules']['kite']['bus_file']}", ''),
            chemistry=config['modules']['kite']['chemistry']

        resources:
            **get_resources("run_kallisto")

        conda: "../envs/kallisto.yaml"

        envmodules: "kallisto/0.46.1"

        shell:
            """
            kallisto bus -i {input[0]} -o {params.output_pref} -x {params.chemistry} -t 8 {input.fastq_files}
            """


if SINGLE_HTO:
    rule run_bustools_correct:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_file']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_file_corrected']}"

        params:
            whitelist=config['references']['whitelist']

        resources:
            **get_resources("run_bustools_correct")

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            bustools correct -w {params.whitelist} {input} -o {output}
            """

elif MULTI_HTO:
    rule run_bustools_correct:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_file']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_file_corrected']}"

        params:
            whitelist=config['references']['whitelist']

        resources:
            **get_resources("run_bustools_correct")

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            bustools correct -w {params.whitelist} {input} -o {output}
            """


if SINGLE_HTO:
    rule run_bustools_sort:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_file_corrected']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_file_sorted']}"

        resources:
            **get_resources("run_bustools_sort")

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            bustools sort -t 4 -o {output} {input}
            """

elif MULTI_HTO:
    rule run_bustools_sort:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_file_corrected']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_file_sorted']}"

        resources:
            **get_resources("run_bustools_sort")

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            bustools sort -t 4 -o {output} {input}
            """


if SINGLE_HTO:
    rule run_bustools_count:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_file_sorted']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['ec_matrix']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['tx']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['features_mismatch_t2g']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_mtx']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_genes']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_barcodes']}"

        params:
            output_pref=lambda wildcards, output: output[0].replace(f"{config['modules']['kite']['bus_count_mtx']}", '')

        resources:
            **get_resources("run_bustools_count")

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            bustools count -o {params.output_pref} --genecounts -g {input[3]} -e {input[1]} -t {input[2]} {input[0]}
            """

elif MULTI_HTO:
    rule run_bustools_count:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_file_sorted']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['ec_matrix']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['tx']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['features_mismatch_t2g']}"

        # priority: 10

        output:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_mtx']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_genes']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_barcodes']}"

        params:
            output_pref=lambda wildcards, output: output[0].replace(f"{config['modules']['kite']['bus_count_mtx']}", '')

        resources:
            **get_resources("run_bustools_count")

        conda: "../envs/bustools.yaml"

        envmodules: "bustools/0.40.0"

        shell:
            """
            bustools count -o {params.output_pref} --genecounts -g {input[3]} -e {input[1]} -t {input[2]} {input[0]}
            """


if SINGLE_HTO:
    rule create_h5ad_bustools:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_mtx']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_genes']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_barcodes']}"

        # priority: 10

        output:
            f"{config['modules']['hashsolo']['h5ad_bustools_dir']}{config['modules']['hashsolo']['fold_struct']}{config['modules']['hashsolo']['bustools_h5ad']}"

        resources:
            **get_resources("create_h5ad_bustools")

        conda: "../envs/basic_sctools.yaml"
        
        shell:
            """
            python3 helper_py_scripts/create_h5ad_from_bustools.py {input[0]} {input[1]} {input[2]} -o {output}
            """

elif MULTI_HTO:
    rule create_h5ad_bustools:
        input:
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_mtx']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_genes']}",
            f"{config['modules']['kite']['kallisto_bustools_dir']}{config['modules']['kite']['fold_struct']}{{i}}{config['modules']['kite']['bus_count_dir']}{config['modules']['kite']['bus_count_barcodes']}"

        # priority: 10

        output:
            f"{config['modules']['hashsolo']['h5ad_bustools_dir']}{config['modules']['hashsolo']['fold_struct']}_{{i}}{config['modules']['hashsolo']['bustools_h5ad']}"

        resources:
            **get_resources("create_h5ad_bustools")

        conda: "../envs/basic_sctools.yaml"
        
        shell:
            """
            python3 helper_py_scripts/create_h5ad_from_bustools.py {input[0]} {input[1]} {input[2]} -o {output}
            """
