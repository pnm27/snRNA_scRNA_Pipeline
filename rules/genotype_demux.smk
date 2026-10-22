# This rule will create multiple runs of cellSNP and vireoSNP for each sample
# Provided a file that contains multiple vcfs per each sample
# present as config['modules']['gt_demux']['vcf_info']
from lib.io import (
    count_matrix,
    previous_solo_h5ad,
    get_cellsnp_inputs,
    get_vir_inputs
)
from lib.utils import get_demux_type
from lib.params import (
    get_params_cics, 
    get_procs_csnp, 
    get_cmd_str_vireo, 
    get_umiTag
)
from functools import partial



# Multi_module branch related function
# def multi_vcfs(wildcards):
#     temp_df = read_files_ext(config['modules']['gt_demux']['vcf_info'])
#     wet_lab_df = read_files_ext(config['experiment']['inputs']['wet_lab_info'])
#     samp = wildcards.id1.replace('-', '_')+'_cDNA'
#     cols1 = config['modules']['gt_demux']['vcf_info_columns']
#     cols2 = config['modules']['gt_demux']['wet_lab_file_cols']
#     set_num = wet_lab_df.loc[wet_lab_df[cols2[1]].str.lower() == samp.lower(), cols2[0]].values[0]
#     col = wildcards.vcf_type

#     return temp_df.loc[temp_df[cols1[0]] == set_num, col].values[0]

# Modify this rule's input to accomodate the basic filtering that 
# is used for calico_solo's input


# Cells come from the hashsolo h5ad that add_vireo adds to, else from the
# count matrix (filtered as in calico_solo's rule)
rule create_inp_cellSNP:
    input:
        branch(
            get_demux_type(config) == 'add_vireo',
            then=partial(previous_solo_h5ad, config=config),
            otherwise=partial(count_matrix, config=config),
        )

    priority: 8

    params:
        extra=partial(get_params_cics, config=config, input=input)

    resources:
        **get_resources("create_inp_cellSNP")

    output:
        f"{config['modules']['gt_demux']['inp_for_cellsnp_dir']}{config['modules']['gt_demux']['fold_struct_filt_bc']}.txt"
    
    conda: "../envs/basic_sctools.yaml"

    shell:
        """
        python3 helper_py_scripts/create_inp_cellSNP.py {input} \
            -o {output} {params.extra}
        sleep 60
        """




# UMI tag is turned on. Therefore, PCR duplicates are included
rule cellSNP:
    input:
        unpack(partial(get_cellsnp_inputs, config=config))

    output:
        f"{config['modules']['gt_demux']['cellsnp_dir']}{config['modules']['gt_demux']['fold_struct']}{config['modules']['gt_demux']['cellsnp_cells']}",
        f"{config['modules']['gt_demux']['cellsnp_dir']}{config['modules']['gt_demux']['fold_struct']}{config['modules']['gt_demux']['cellsnp_base']}"

    params:
        # ref_snps=config['modules']['gt_demux']['ref_snps'],
        umi_tag=partial(get_umiTag, config=config),
        cell_tag=config['modules']['gt_demux']['cell_tag'],
        processors=partial(get_procs_csnp, config=config),
        min_maf=config['modules']['gt_demux']['min_maf'],
        min_ct=config['modules']['gt_demux']['min_aggr_count'],
        output_prefix=lambda _, output: (
            output[0].replace(
                f"/{config['modules']['gt_demux']['cellsnp_cells']}", 
                ''
            )
        ),
        threads=config['modules']['gt_demux']['bcftools_thread']


    resources:
        **get_resources("cellSNP")
    
    conda: "../envs/gt_demux.yaml"

    envmodules: "bcftools/1.15.1"

    shell:
        """
        set -x
        cellsnp-lite -b {input.barcodesFile} -s {input.bam} \
            -R {input.regionsFile} -O {params.output_prefix} \
            -p {params.processors} --minMAF {params.min_maf} \
            --minCOUNT {params.min_ct} --cellTAG {params.cell_tag} \
            --UMItag {params.umi_tag} --genotype --gzip
        set +x
        """


rule vireoSNP:
    input:
        unpack(partial(get_vir_inputs, config=config))
        
    output:
        (
            f"{config['modules']['gt_demux']['vireosnp_dir']}"
            f"{config['modules']['gt_demux']['fold_struct']}"
            f"{config['modules']['gt_demux']['donors_classification']}"
        )

    # Read by demul_samples.py (next to donor_ids.tsv) for variant matching stats
    log:
        (
            f"{config['modules']['gt_demux']['vireosnp_dir']}"
            f"{config['modules']['gt_demux']['fold_struct']}"
            "vireo.log"
        )

    params:
        geno_tag=config['modules']['gt_demux']['donor_genotype'],
        output_prefix=lambda _, output: (
            output[0].replace(
                f"/{config['modules']['gt_demux']['donors_classification']}", 
                ''
            )
        ),
        cmd_str=partial(get_cmd_str_vireo, config=config)

    resources:
        **get_resources("vireoSNP")
        
    conda: "../envs/gt_demux.yaml"

    shell:
        """        
        set -x
        vireo -c {input.cellsnpCells} -o {params.output_prefix} \
            -t {params.geno_tag} --noPlot --randSeed 100 \
            {params.cmd_str} 2>&1 | tee {log}
        set +x
        """
