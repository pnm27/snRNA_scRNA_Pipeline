# This rule will create multiple runs of cellSNP and vireoSNP for each sample
# Provided a file that contains multiple vcfs per each sample
# present as config['gt_demux_pipeline']['vcf_info']
from lib.utils import read_files_ext, ret_cols
from lib.io import (
    get_filt_barcodes, 
    get_cellsnp_inputs,
    get_vir_inputs
)
from lib.params import (
    get_params, 
    get_cmd_str_csnp, 
    get_cmd_str_vireo, 
    get_umiTag
)



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
    return 1000+500*(attempt-1)


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
        lambda wc: get_filt_barcodes(wc, config, global_vars)

    priority: 8

    params:
        # DEPRACATED
        # inp_type=get_inp_type, 
        # When a run of calico_solo exists
        # col_name=config['gt_demux_pipeline']['demux_col'], # Name of the anndata's obs column that contains classification of cells
        # bc_len=config['gt_demux_pipeline']['barcode_len'], # Barcode length
        # keep_all_cells=config['gt_demux_pipeline']['include_all_cells'], # Include all cells (don't remove cells prev classified as doublets, etc.)
        # doub=config['gt_demux_pipeline']['doublet'], # Doublets named as
        # neg=config['gt_demux_pipeline']['negative'], # Negatives named as
        # na=config['gt_demux_pipeline']['na'], # Cells not present in hashsolo named as
        # When no previous runs of calico_solo exists
        # mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        # min_genes=config['min_genes_per_cell'], # Min #genes per cell
        # min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        # genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
        # mito_prefix=config['mito_prefix'], # Mitochondrial genes' (names') prefix
        extra=lambda wc, input: get_params(wc, input, config, global_vars)

    resources:
        cpus_per_task=2, # For snakemake > v8
        mem_mb=allocate_mem_CICS,
        time_min=allocate_time_CICS

    # group: "genotype-demux"
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
        lambda wc: get_cellsnp_inputs(wc, config)

    # group: "genotype-demux"

    output:
        f"{config['gt_demux_pipeline']['cellsnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['cellsnp_cells']}",
        f"{config['gt_demux_pipeline']['cellsnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['cellsnp_base']}"

    params:
        # ref_snps=config['gt_demux_pipeline']['ref_snps'],
        umi_tag=lambda wc: get_umiTag(wc, config),
        cell_tag=config['gt_demux_pipeline']['cell_tag'],
        processors=config['gt_demux_pipeline']['n_proc'],
        min_maf=config['gt_demux_pipeline']['min_maf'],
        min_ct=config['gt_demux_pipeline']['min_aggr_count'],
        output_prefix=lambda wildcards, output: output[0].replace(f"/{config['gt_demux_pipeline']['cellsnp_cells']}", ''),
        filt_vcf_dir=f"{config['gt_demux_pipeline']['filt_vcf_dir']}{config['fold_struct_gt_demux']}"[:-1], # remove trailing forward slash
        threads=config['gt_demux_pipeline']['bcftools_thread'],
        cmd_str=lambda wc, input: get_cmd_str_csnp(wc, input)

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
        read -r -a array <<< "{input}"
        n_procs=$(( {params.processors} >  ( {resources.cpus_per_task} * 2 ) ? {params.processors} : ( {resources.cpus_per_task} * 2 ) ))
        if [[ "${{#array[@]}}" -lt 4 ]]; then
            set -x
            cellsnp-lite {params.cmd_str} -R ${{array[2]}} -O {params.output_prefix} \
                -p {params.processors} --minMAF {params.min_maf} \
                --minCOUNT {params.min_ct} --cellTAG {params.cell_tag} \
                --UMItag {params.umi_tag} --genotype --gzip
        else
            set -x
            if [ ! -f {params.filt_vcf_dir}"/0002.vcf.gz" ]; then
                bcftools isec --threads {params.threads} -e- -i'INFO/AF>0.25' \
                    -Oz -p {params.filt_vcf_dir} ${{array[@]: -2:2}}
            fi
            cellsnp-lite {params.cmd_str} -O {params.output_prefix} \
                -R {params.filt_vcf_dir}"/0002.vcf.gz" -p {params.processors} \
                --minMAF {params.min_maf} --minCOUNT {params.min_ct} \
                --cellTAG {params.cell_tag} --UMItag {params.umi_tag} \
                --genotype --gzip
        fi
        set +x
        """


rule vireoSNP:
    input:
        lambda wc: get_vir_inputs(wc, config)
        
    output:
        f"{config['gt_demux_pipeline']['vireosnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['donors_classification']}"

    # group: "genotype-demux"

    params:
        # donor_info=get_donor_info,
        geno_tag=config['gt_demux_pipeline']['donor_genotype'],
        output_prefix=lambda wildcards, output: output[0].replace(f"/{config['gt_demux_pipeline']['donors_classification']}", ''),
        cmd_str=lambda wc, input: get_cmd_str_vireo(wc, input, config)

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
        vireo {params.cmd_str} -o {params.output_prefix} -t {params.geno_tag} \
            --noPlot --randSeed 100
        set +x
        """
