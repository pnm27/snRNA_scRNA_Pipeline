def get_bam_inputs(wildcards):
    
    if config['identify_swaps']['mbv_inp'] == 'vireo_outs':
        if 'multiome' in config['last_step'].lower():
            return (
                f"{config['split_bams_pipeline']['split_bams_dir2']}"
                f"{config['fold_struct_bam_split2']}"
                f"{config['fold_struct_gt_demux_redo']}.bam"
            ).format(**wildcards)
        else:
            return (
                f"{config['split_bams_pipeline']['split_bams_dir2']}"
                f"{config['fold_struct_bam_split2']}"
                f"{config['fold_struct_gt_demux_redo']}.bam"
            )
    else:
        return config['identify_swaps']['mbv_inp']


# Resource Allocation ------------------
def allocate_mem_QM(wildcards):
    return 5000


def allocate_time_QM(wildcards, attempt):
    return 60 + 20*(attempt-1)

# --------------------------------------

rule qtltools_mbv:
    input:
        bam=get_bam_inputs,
        ref_snps=config['identify_swaps']['ref_vcf']

    output:
        f"{config['identify_swaps']['mbv_out_dir']}{config['fold_struct_swaps_check']}{config['identify_swaps']['mbv_suffix']}"
        

    params:
        extra_params=config['identify_swaps']['mbv_extra_opt']

    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_QM,
        time_min=allocate_time_QM

    envmodules:
        "qtltools/1.3"

    # Use one of the below
    wildcard_constraints:
        # donor=r"(?<=cDNA_|ATAC_).+" #WILDCARDS # For multiome
        donor=r"(?:.*)(?<=\/)([^/]*)" #WILDCARDS # Match everything except the last '/'

    shell:
        """
        if [[ "{params.extra_params}" == "None" ]]; then
            QTLtools mbv --bam {input.bam} --out {output} --vcf {input.ref_snps}
        else
            QTLtools mbv --bam {input.bam} --out {output} --vcf {input.ref_snps} {params.extra_params}
        fi
        
        """