def get_bam_inputs(wildcards):
    if config['identify_swaps']['mbv_inp'] == 'vireo_outs':
        return f"{config['split_bams_pipeline_gt_demux']['split_bams_dir2']}{config['fold_struct_bam_split2']}{config['fold_struct_gt_demux_redo']}.bam"
    else:
        return config['identify_swaps']['mbv_inp']



rule qtltools_mbv:
    input:
        bam=get_bam_inputs,
        ref_snps=config['identify_swaps']['ref_vcf']

    output:
        f"{config['identify_swaps']['mbv_out_dir']}{config['fold_struct_swaps_check']}{config['identify_swaps']['mbv_suffix']}"
        

    params:
        extra_params=config['identify_swaps']['mbv_extra_opt']

    threads: 1

    resources:
        mem_mb=allocate_mem_QM,
        time_min=allocate_time_QM

    envmodules:
        "qtltools/1.3"

    shell:
        """
        ml qtltools/1.3
        if [[ "{params.extra_params}" == "None" ]]; then
            QTLtools mbv --bam {input.bam} --out {output} --vcf {input.ref_snps}
        else
            QTLtools mbv --bam {input.bam} --out {output} --vcf {input.ref_snps} {params.extra_params}
        fi
        
        """