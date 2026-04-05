import pandas as pd, os
from typing import Literal


# Run STARsolo_sort, index_bams and filter_GeneFull_STARsolo
def targets_STARsolo() -> str:

    inp_pref = config['STARsolo_pipeline']['bams_dir']
    bai_suff = config['STARsolo_pipeline']['bai']
    # gf_mat = config['STARsolo_pipeline']['genefull_matrix']
    # gf_features = config['STARsolo_pipeline']['genefull_features']
    # barcodes_stats = config['STARsolo_pipeline']['barcodes_stats']
    f_st = config['fold_struct']

    # target_list = [f"{inp_pref}{f_st}{bai_suff}", f"{inp_pref}{f_st}{gf_mat}", 
    #     f"{inp_pref}{f_st}{gf_features}", f"{inp_pref}{f_st}{barcodes_stats}"]
    # return target_list
    return f"{inp_pref}{f_st}{bai_suff}"


# Add support for multiome
# Add CollectInsertSizeMetrics for ATAC part of multiome
def targets_PICARD(
    progs: Literal['all', 'gc', 'rnaseq'] = 'all'
    ) -> str:
    
    inp_pref = config['STARsolo_pipeline']['bams_dir']
    rna_seq_suff = config['picard_pipeline']['rnaseq_metrics']
    # gc_met_suff = config['picard_pipeline']['gc_bias_metrics']
    gc_summ_suff = config['picard_pipeline']['gc_summary_metrics']
    f_st=config['fold_struct']

    files_dict = {
        'rnaseq':rna_seq_suff, 
        # 'gc':[gc_met_suff, gc_summ_suff],
        'gc': gc_summ_suff
        }
    target_list = []
    if progs == 'all':
       # Run both PICARD programs
        # for k, val in files_dict.items():
            # if isinstance(val, list):
            #     for l_val in val:
            #         target_list.append(f"{inp_pref}{f_st}{l_val}")

            # else:
        target_list.extend([f"{inp_pref}{f_st}{val}" for k, val in files_dict.items() ])

    else:     
        # target_list.extend( [f"{inp_pref}{f_st}{val}"  for val in files_dict[progs]] )
        val = files_dict[progs]
        target_list.append(f"{inp_pref}{f_st}{val}")

    # elif progs == "gc":
    #     target_list.extend( [f"{inp_pref}{f_st}{val}"  for val in files_dict[progs]] )
    # else:
    #     return []

    return ','.join(target_list)



# Temporary function, will change to targets_gt_demux
def targets_cellSNP(multiome: bool = False) -> str:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    out_dir = config['gt_demux_pipeline']['cellsnp_dir']
    fs_gt = config['fold_struct_gt_demux']
    cells_suff = config['gt_demux_pipeline']['cellsnp_cells']
    # base_suff = config['gt_demux_pipeline']['cellsnp_base']
    for d in sub_dir:
        # target_list.extend([
        #     os.path.join(f"{out_dir}{fs_gt}", d, f"{cells_suff}"), 
        #     os.path.join(f"{out_dir}{fs_gt}", d, f"{base_suff}")
        #     ])
        target_list.append(os.path.join(f"{out_dir}{fs_gt}", d, f"{cells_suff}"))

    # STARsolo* + PICARD (any) progs
    # if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
        # target_list.extend(targets_PICARD( progs=progs))
        # return ','.join(target_list) + ',' + targets_PICARD(progs=progs)

    return ','.join(target_list)



def targets_gt_demux(h5ad: bool = False, 
    multiome: bool = False) -> str:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    if not h5ad:
        for d in sub_dir:
        
            out_dir = config['gt_demux_pipeline']['vireosnp_dir']
            fs_gt = config['fold_struct_gt_demux']
            suff = config['gt_demux_pipeline']['donors_classification']
            target_list.append(os.path.join(f"{out_dir}{fs_gt}", d ,f"{suff}"))
    # For multiome combine cDNA and ATAC based demultiplexing into
    # a single h5ad
    else:
        out_dir = config['gt_demux_pipeline']['final_count_matrix_dir']
        fs_gt = config['fold_struct_demux']
        suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']
        target_list.append(os.path.join(f"{out_dir}{fs_gt}{suff}"))


    # STARsolo* + PICARD (any) progs
    # if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
    #     # target_list.extend(targets_PICARD( progs=progs))
    #     return ','.join(target_list) + ',' + targets_PICARD(progs=progs)

    return ','.join(target_list)


# For multi-vcf_inputs
# CONSIDER RE-EVALUATE
def targets_gt_demux2(multiome: bool = False) -> str:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    for d in sub_dir:
        out_dir = config['gt_demux_pipeline']['final_count_matrix_dir']
        fs_gt = config['fold_struct_demux']
        
        # target_list = [f"{out_dir}{fs_gt}" ]
        target_list.append(os.path.join(f"{out_dir}", d, f"{fs_gt}"))
    # STARsolo* + PICARD (any) progs
    # if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
    #     return ','.join(target_list) + ',' + targets_PICARD(progs=progs)

    # return ','.join(target_list)


# NEED TO WORK
def targets_SplitBams(multiome: bool = False) -> str:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    for d in sub_dir:
        if config['gt_check']:
            out_dir = config['split_bams_pipeline']['split_bams_proxy_dir2']
        else:
            out_dir = config['split_bams_pipeline']['split_bams_proxy_dir']
        op = config['fold_struct_bam_split1']

        f = os.path.join(f"{out_dir}", f"{op}") if d == '' \
            else os.path.join(f"{out_dir}", f"{op}", d)
        target_list.append(f)
    
    # STARsolo* + PICARD (any) progs
    # if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
    #     target_list.extend(targets_PICARD( progs=progs, 
    #         ))
    
    # return target_list


def targets_resolve_swaps_gt_demux() -> str:
    out_dir = config['finalize_demux']['out_dir']
    fs = config['fold_struct_gt_demux_redo']
    dons_file = config['finalize_demux']['donor_info_swap_file']
    suff = config['gt_demux_pipeline']['donors_classification']
    
    # return [f"{out_dir}{fs}{suff}"]
    return f"{out_dir}{fs}{suff}"


def targets_gt_demux_identify_swaps(multiome: bool = False) -> str:

    out_dir = config['identify_swaps']['mbv_out_dir']
    fs = config['fold_struct_swaps_check']
    suff = config['identify_swaps']['mbv_suffix']
    # target_list = []  
    # target_list.append(os.path.join(f"{out_dir}", f"{fs}{suff}"))

    return os.path.join(f"{out_dir}", f"{fs}{suff}")


def targets_multibamsummary() -> list:
    out_dir = config['multiBamSummary']['outdir']
    fs = config['fold_struct_deeptools']
    suff = ".npz"
    
    return [f"{out_dir}{fs}{suff}"]


def targets_multibamsummaryPlotCorr() -> list:
    out_dir = config['plotCorrelation']['outdir']
    fs = config['fold_struct_deeptools']
    suff = config['plotCorrelation']['out_fmt']
    
    return [f"{out_dir}{fs}{suff}"]


# Need to add support for PICARD metrics
# and create h5ad
def targets_multiome(last, h5ad=False, multiome=False) -> list:
    if last == 'alignment':
        out_dir = config['cellranger_arc_count']['bams_dir']
        return [
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/barcodes.tsv.gz",
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/features.tsv.gz",
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/matrix.mtx.gz",
        ]
    elif last == 'vireo':
        return targets_gt_demux(h5ad=h5ad, multiome=multiome)

    elif last == 'splitBams':
        return targets_SplitBams(multiome=multiome)

    elif last == 'identifySwaps':
        return targets_gt_demux_identify_swaps(multiome=multiome)


# To run STARsolo* + kb pipeline + (optional)PICARD progs
# MIGHT BE REMOVED
def targets_all(
    progs: Literal['all', 'gc', 'rnaseq', None] = None
    ) -> list:
    
    demuxed_mat_dir = config['hashsolo_demux_pipeline']['final_count_matrix_dir']
    demuxed_info_dir = config['hashsolo_demux_pipeline']['demultiplex_info_dir']    
    demux_mat_suff = config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']
    demux_info_suff = config['hashsolo_demux_pipeline']['demultiplex_info']
    f_st = config['fold_struct_demux']

    # STARsolo* + PICARD (any) progs
    if progs is not None:
        target_list = targets_PICARD(progs=progs)

    else: # For any wrong values to progs or progs == None or just PICARD == False or everything else, run STARsolo*
        target_list = targets_STARsolo()


    # kb pipeline
    target_list.extend( [f"{demuxed_mat_dir}{f_st}{demux_mat_suff}", f"{demuxed_info_dir}{f_st}{demux_info_suff}"] )
    
    return target_list


# Add wildcards info here in the 'expand' function
def produce_targets(wildcards) -> list:
    target_step = config['last_step'].lower()
    
    if target_step is None:
        raise ValueError("The variable last_step in config file is set to None!")

    metrics = config['picard_metrics'].lower() if config['picard_metrics'] is not None else None
    create_h5ad = config['create_final_h5ad'] if 'gt_demux' in target_step else False
    multiome = True if 'multiome' in target_step else False

    final_target_list = []
    # for pool in POOLS:
    #     wc = {"pool": pool}
        # DEPRACATED
        # if target_step == "all":
        #     target_files = targets_all(progs=metrics)
        #     target_files.extend(targets_gt_demux(progs=metrics))
        #     final_target_list = [expand(f"{target}", zip, **wc_d) for target in target_files] # if multiple wildcards are present then either 'zip' or 'product' (default)
        # NEEDS RE-CONSIDERATION
    if target_step == "all_multi_vcf":
        target_files = targets_all(progs=metrics)
        final_target_list = [expand(f"{target}", zip, **wc_d) for target in target_files] # if multiple wildcards are present then either 'zip' or 'product' (default)
        # Adding gt_demux_multi_vcf_options
        target_files = targets_gt_demux2()
        suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']

        # If multiple vcfs per sample needs to be run
        if isinstance(VCF_TYPE, list):
            
            temp_list=[expand(f"{target}", zip, **wc_d) for target in target_files][0] # Multiple wildcards example
            for id, target in enumerate(target_files):
                if id == 0:
                    temp_list = expand(f"{target}", zip, **wc_d)
                    final_target_list.extend(
                        [expand(f"{target}_{{vcf_type}}{suff}", vcf_type=VCF_TYPE) for t in temp_list]
                        )
                else:
                    final_target_list.extend(expand(f"{target}", zip, **wc_d))

        # Single vcf out of a set of multiple vcfs
        elif isinstance(VCF_TYPE, str):
            
            final_target_list.extend(
                [expand(f"{target}_{VCF_TYPE}{suff}", zip, **wc_d) for target in target_files]
                ) # Multiple wildcards example
            # final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", **wc_d) for target in target_files] # Single wildcard
        else:
            raise ValueError("The wildcard used to test multi_vcf "
                "input (VCF_TYPE) is of unexpected type! Please check")

    # NEED to RE-EVALUATE
    elif target_step == "starsolo_kb_solo":
        target_files = targets_all(progs=metrics)
        final_target_list = [expand(f"{target}", zip, **wc_d) for target in target_files]     

    elif target_step == "starsolo_cellsnp":
        target_files = targets_cellSNP(multiome=False)
        final_target_list.append(target_files.format(pool=pool))
        if metrics is not None:
            final_target_list += [
                target.format(pool=pool) 
                for target in targets_PICARD(progs=metrics).split(',')
                ]
        
    elif target_step == "starsolo_gt_demux":
        target_files = targets_gt_demux(h5ad=create_h5ad, 
                    multiome=multiome)
        final_target_list.append(target_files.format(pool=pool))
        if metrics is not None:
            final_target_list += [
                target.format(pool=pool) 
                for target in targets_PICARD(progs=metrics).split(',')
                ]


    # Multi_vcf inputs or when a subdir with the vcf name is needed
    # NEEDS RE-CONSIDERATION
    elif target_step == "starsolo_gt_demux_multi_vcf":
        target_files = targets_gt_demux2()
        suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']

        # If multiple vcfs per sample needs to be run
        if isinstance(VCF_TYPE, list):
            
            temp_list= [expand(f"{target}", zip, **wc_d) for target in target_files][0] # Multiple wildcards example
            for id, target in enumerate(target_files):
                if id == 0:
                    temp_list = expand(f"{target}", zip, **wc_d)
                    final_target_list.extend(
                        [expand(f"{target}_{{vcf_type}}{suff}", zip, **wc_d, vcf_type=VCF_TYPE) for t in temp_list]
                        )
                else:
                    final_target_list.extend(expand(f"{target}", zip, **wc_d))

            # temp_list= [expand(f"{target}", **wc_d) for target in target_files][0] # Single wildcard
            # final_target_list=[expand(f"{t}_{{vcf_type}}{suff}", vcf_type=VCF_TYPE) for t in temp_list]

        # Single vcf out of a set of multiple vcfs
        elif isinstance(VCF_TYPE, str):
            
            final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", zip, **wc_d) for target in target_files] # Multiple wildcards example
            # final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", **wc_d) for target in target_files] # Single wildcard
        else:
            raise ValueError("The wildcard used to test multi_vcf input (VCF_TYPE) is of unexpected type! Please check")

    elif target_step == "multibamsummaryplotcorr":
        # global_vars.ONLY_SOLO = True
        target_files = targets_multibamsummaryPlotCorr()
        final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files[0]]

    elif target_step == "multibamsummary":
        # global_vars.ONLY_SOLO = True
        target_files = targets_multibamsummary()
        final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files[0]]

    elif target_step == "multiome_alignment":
        target_files = targets_multiome(last="alignment")
        final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

    # Support creation of h5ad or not
    elif target_step == "multiome_gt_demux":
        target_files = targets_multiome(last="vireo", progs=metrics, 
        multiome=multiome, h5ad=create_h5ad)
        final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

    elif target_step == "multiome_split_bams_gt_demux":
        target_files = targets_multiome(last="splitBams", progs=metrics,
        multiome=multiome, h5ad=create_h5ad)
        suff = ".txt"
        
        final_target_list = []

        for id, target in enumerate(target_files):
            if id <= 1:
                final_target_list.extend(expand(f"{target}{suff}", zip, **wc_d))
            else:
                final_target_list.extend(expand(f"{target}", zip, **wc_d))

    elif target_step == "multiome_gt_demux_identify_swaps":
        target_files = targets_multiome(last="identifySwaps", progs=metrics,
        multiome=multiome, h5ad=create_h5ad)
        final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

    # NEED to RE-EVALUATE
    elif "starsolo" in target_step:
        if metrics is None:
            target_files = targets_STARsolo()
            final_target_list = [target_files.format(pool=pool)]
        else:
            target_files = targets_PICARD(progs=metrics)
            final_target_list = target_files.split(',')
            final_target_list = [f"{target}".format(pool=pool) for target in final_target_list]

    elif any([ 
        f in config['last_step'].lower() 
        for f in ['split_bams', 'identify_swaps']
        ]):
        expand("resolved/{pool}.done", pool=POOLS)

    return final_target_list


# Add wildcards info here in the 'expand' function
def produce_targets_dynamic(wildcards) -> list:
    target_step = config['last_step'].lower()
    final_target_list = []
    # for pool in POOLS:
    #     wc = {"pool": pool}
    if target_step == "starsolo_split_bams":
        target_files = targets_SplitBams()
        # final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]
        wc_dict = get_real_donors({"pool": wildcards.pool})
        return expand(target_files, zip, **wc_dict)

    # RE-EVALUTE
    elif target_step == "starsolo_split_bams_gt_demux":
        target_files = targets_SplitBams()
        suff = ".txt"
        # final_target_list= [expand(f"{target}{suff}", zip, **wc_d) for target in target_files][0] # Single wildcard
        wc_dict = get_real_donors({"pool": wildcards.pool})
        final_target_list = []
        for id, target in enumerate(target_files):
            if id == 0:
                final_target_list.extend(expand(f"{target}{suff}", zip, **wc_dict))
            else:
                final_target_list.extend(expand(f"{target}", zip, **wc_dict))

    # RE-EVALUTE
    # Multi_vcf inputs or when a subdir with the vcf name is needed
    elif target_step == "starsolo_split_bams_gt_demux_multi_vcf":
        # Add final steps for split_bams
        target_files = targets_SplitBams()
        suff = ".txt"
        # If multiple vcfs per sample needs to be run or not
        if isinstance(VCF_TYPE, list):
            
            temp_list = [expand(f"{target}", zip, **wc_d) for target in target_files][0] # Multiple wildcards example
            # temp_list= [expand(f"{target}", **wc_d) for target in target_files][0] # Single wildcard
            return [expand(f"{t}_{{vcf_type}}{suff}", vcf_type=VCF_TYPE) for t in temp_list]

        # Single vcf out of a set of multiple vcfs
        elif isinstance(VCF_TYPE, str):
            
            return [expand(f"{target}_{VCF_TYPE}{suff}", zip, **wc_d) for target in target_files] # Multiple wildcards example
            # final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", **wc_d) for target in target_files] # Single wildcard
        else:
            raise ValueError("The wildcard used to test multi_vcf input (VCF_TYPE) is of unexpected type! Please check")

        # If multiple vcfs per sample needs to be run
        if isinstance(VCF_TYPE, list):
            
            # temp_list = [expand(f"{target}", zip, **wc_d) for target in target_files][0] # Multiple wildcards example
            for id, target in enumerate(target_files):
                if id == 0:
                    temp_list = [expand(f"{target}", zip, **wc_d)]
                    final_target_list.extend(
                        [expand(os.path.join(f"{t}", f"{{vcf_type}}{suff}"), vcf_type=VCF_TYPE) \
                        for t in temp_list]
                    )
                else:
                    final_target_list.extend(expand(f"{target}", zip, **wc_d))
            # temp_list= [expand(f"{target}", **wc_d) for target in target_files][0] # Single wildcard
            

        # Single vcf out of a set of multiple vcfs
        elif isinstance(VCF_TYPE, str):

            for id, target in enumerate(target_files):
                if id == 0:
                    temp_list = expand(f"{target}", zip, **wc_d)
                    final_target_list.extend([
                        os.path.join(f"{target}", f"{VCF_TYPE}{suff}") \
                        for t in temp_list])
                else:
                    final_target_list.extend(expand(f"{target}", zip, **wc_d))

            # final_target_list2 = [expand(f"{target}_{VCF_TYPE}{suff}", zip, **wc_d) for target in target_files] # Multiple wildcards example
            # final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", **wc_d) for target in target_files] # Single wildcard
        else:
            raise ValueError("The wildcard used to test multi_vcf input (VCF_TYPE) is of unexpected type! Please check")


    elif target_step == "starsolo_gt_demux_identify_swaps":
        target_files = targets_gt_demux_identify_swaps(multiome=False)
        wc_dict = get_real_donors({"pool": wildcards.pool})
        return  expand(target_files, zip, **wc_dict)
        # final_target_list += expand(target_files, zip, **wc_dict)
        # if metrics is not None:
        #     final_target_list += [
        #         target.format(pool=pool) 
        #         for target in targets_PICARD(progs=metrics).split(',')
        #         ]
    
    # MIGHT REMOVE
    elif target_step == "starsolo_resolve_swaps_gt_demux":
        target_files = targets_resolve_swaps_gt_demux()
        return [ expand(f"{target}", zip, **wc_d) for target in target_files]

    # return final_target_list

