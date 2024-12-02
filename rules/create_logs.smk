# Add wildcards info here in the 'expand' function
def stats_produce_inp(wildcards):
    STAR_log=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['STAR_log_final']}", id1=sample_name)
    SS_G_Feat=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['gene_features']}", id1=sample_name)
    SS_GF_Feat=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_features']}", id1=sample_name)
    SS_G_Summ=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['gene_summary']}", id1=sample_name)
    SS_GF_Summ=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_summary']}", id1=sample_name)
    SS_Barcodes=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['barcodes_stats']}", id1=sample_name)
    PICARD_GC=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['gc_summary_metrics']}", id1=sample_name)
    PICARD_RNAseq=expand(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['rnaseq_metrics']}", id1=sample_name)
    Demultiplex_info=expand(f"{config['hashsolo_demux_pipeline']['demultiplex_info_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['demultiplex_info']}", id1=sample_name)

    if config['last_step'] == "all":
        
        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_GC, PICARD_RNAseq, Demultiplex_info]

    elif config['last_step'] == "STARsolo_kb_solo":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, Demultiplex_info]

    elif config['last_step'] == "STARsolo":
        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes]

    elif config['last_step'] == "STARsolo_PICARD":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_GC, PICARD_RNAseq]

    elif config['last_step'] == "STARsolo_rnaseqmet":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_RNAseq]

    elif config['last_step'] == "STARsolo_gcbiasmet":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_GC]

    # Not yet implemented
    elif config['last_step'] == "STARsolo_gt_demux":
        pass

    else:
        raise ValueError(f"Invalid Value to {config['last_step']} in the yaml file")



def stats_produce_params(wildcards, input):

    # Files for update log files
    log_dict = {'STAR_log':[config['STARsolo_pipeline']['STAR_log_final'], '--ss_l'], 'PICARD_GC':[config['picard_pipeline']['gc_summary_metrics'], '--pc_gc'],
                'PICARD_RNAseq':[config['picard_pipeline']['rnaseq_metrics'], '--pc_rs'],
                'SS_G_Feat':[config['STARsolo_pipeline']['gene_features'], '--ss_g_f'], 'SS_GF_Feat':[config['STARsolo_pipeline']['genefull_features'], '--ss_gf_f'],
                'SS_G_Summ':[config['STARsolo_pipeline']['gene_summary'], '--ss_g_s'], 'SS_GF_Summ':[config['STARsolo_pipeline']['genefull_summary'], '--ss_gf_s'],
                'SS_Barcodes':[config['STARsolo_pipeline']['barcodes_stats'], '--ss_bc'], 'Demultiplex_info':[config['hashsolo_demux_pipeline']['demultiplex_info'], '--dem_info']}

    cons_param = ""
    for i in range(len(input)):
        s_temp = pd.Series(input[i])
        for k, v in log_dict.items():
            # All files should contain the pattern string
            if all(s_temp.str.contains(v[0])):
                curr_param = v[1] + f" {{input[{i}]}}"
                cons_param += curr_param


    return cons_param


# Rule to run update.py
rule update_stats:
    input:
        inp_files = stats_produce_inp
            
    params:
        map_file=config['meta_data'],
        cl_inp_params=stats_produce_params

    output:
        config['log_all_stats']

    threads: 1

    resources:
        mem_mb=allocate_mem_US,
        time_min=allocate_time_US

    shell:
        """
        python3 helper_py_scripts/update_logs.py {params.cl_inp_params} -m {params.map_file} -o {output}
        sleep 60
        """
