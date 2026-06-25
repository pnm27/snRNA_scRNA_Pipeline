# This rule will create multiple runs of cellSNP and vireoSNP for each sample
# Provided a file that contains multiple vcfs per each sample
# present as config['gt_demux_pipeline']['vcf_info']
from lib.io import (
    get_filt_barcodes, 
    get_cellsnp_inputs,
    get_vir_inputs
)
from lib.params import (
    get_params_cics, 
    get_procs_csnp, 
    get_cmd_str_vireo, 
    get_umiTag
)
from functools import partial



# Multi_module branch related function
# def multi_vcfs(wildcards):
#     temp_df = read_files_ext(config['gt_demux_pipeline']['vcf_info'])
#     wet_lab_df = read_files_ext(config['wet_lab_info'])
#     samp = wildcards.id1.replace('-', '_')+'_cDNA'
#     cols1 = config['gt_demux_pipeline']['vcf_info_columns']
#     cols2 = config['gt_demux_pipeline']['wet_lab_file_cols']
#     set_num = wet_lab_df.loc[wet_lab_df[cols2[1]].str.lower() == samp.lower(), cols2[0]].values[0]
#     col = wildcards.vcf_type

#     return temp_df.loc[temp_df[cols1[0]] == set_num, col].values[0]

# Modify this rule's input to accomodate the basic filtering that 
# is used for calico_solo's input


# Resource Allocation ------------------
def allocate_mem_CICS(wildcards, attempt):
    return 2000+500*(attempt-1)


def allocate_time_CICS(wildcards, attempt):
    return 2*attempt+1


def allocate_mem_cS(wildcards, attempt): # WILDCARDS
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 200+80*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 800+100*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 1300+200*(attempt-1)
        else:
            return 1500+200*(attempt-1)
    else:
        return 200+200*(attempt-1)


def allocate_time_cS(wildcards, attempt): # WILDCARDS
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 30+20*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 70+20*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 180+40*(attempt-1)
        else:
            return 240+60*(attempt-1)
    else:
        return 240+60*(attempt-1)


def allocate_mem_vS(wildcards, attempt): # WILDCARDS
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 3000+100*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 5000+400*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 10000+600*(attempt-1)
        else:
            return 22000+500*(attempt-1)
    else:
        return 35000+500*(attempt-1)


def allocate_time_vS(wildcards, attempt): # WILDCARDS
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 10+5*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 70+20*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 180+40*(attempt-1)
        else:
            return 180+60*(attempt-1)
    else:
        return 180+60*(attempt-1)

# --------------------------------------

rule create_inp_cellSNP:
    input:
        partial(get_filt_barcodes, config=config)

    priority: 8

    params:
        extra=partial(get_params_cics, config=config, input=input)

    resources:
        cpus_per_task=2, # For snakemake > v8
        mem_mb=allocate_mem_CICS,
        time_min=allocate_time_CICS

    # For snakemake < v8
    # threads: 2

    output:
        f"{config['gt_demux_pipeline']['inp_for_cellsnp_dir']}{config['fold_struct_filt_bc']}.txt"
    
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
        f"{config['gt_demux_pipeline']['cellsnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['cellsnp_cells']}",
        f"{config['gt_demux_pipeline']['cellsnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['cellsnp_base']}"

    params:
        # ref_snps=config['gt_demux_pipeline']['ref_snps'],
        umi_tag=partial(get_umiTag, config=config),
        cell_tag=config['gt_demux_pipeline']['cell_tag'],
        processors=partial(get_procs_csnp, config=config),
        min_maf=config['gt_demux_pipeline']['min_maf'],
        min_ct=config['gt_demux_pipeline']['min_aggr_count'],
        output_prefix=lambda _, output: (
            output[0].replace(
                f"/{config['gt_demux_pipeline']['cellsnp_cells']}", 
                ''
            )
        ),
        threads=config['gt_demux_pipeline']['bcftools_thread']


    # For snakemake < v8
    # threads: 8

    resources:
        cpus_per_task=8, # For snakemake > v8
        mem_mb=allocate_mem_cS,
        time_min=allocate_time_cS
    
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
            f"{config['gt_demux_pipeline']['vireosnp_dir']}"
            f"{config['fold_struct_gt_demux']}"
            f"{config['gt_demux_pipeline']['donors_classification']}"
        )

    params:
        geno_tag=config['gt_demux_pipeline']['donor_genotype'],
        output_prefix=lambda _, output: (
            output[0].replace(
                f"/{config['gt_demux_pipeline']['donors_classification']}", 
                ''
            )
        ),
        cmd_str=partial(get_cmd_str_vireo, config=config)

    # For snakemake < v8
    # threads: 7

    resources:
        cpus_per_task=7, # For snakemake > v8
        mem_mb=allocate_mem_vS,
        time_min=allocate_time_vS
        
    conda: "../envs/gt_demux.yaml"

    shell:
        """        
        set -x
        vireo -c {input.cellsnpCells} -o {params.output_prefix} \
            -t {params.geno_tag} --noPlot --randSeed 100 \
            {params.cmd_str}
        set +x
        """
