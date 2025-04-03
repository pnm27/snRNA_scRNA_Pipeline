import pandas as pd, os


# Run STARsolo_sort, index_bams and filter_GeneFull_STARsolo
def targets_STARsolo(conf_f) -> list:

    inp_pref = conf_f['STARsolo_pipeline']['bams_dir']
    bai_suff = conf_f['STARsolo_pipeline']['bai']
    gf_mat = conf_f['STARsolo_pipeline']['genefull_matrix']
    gf_features = conf_f['STARsolo_pipeline']['genefull_features']
    barcodes_stats = conf_f['STARsolo_pipeline']['barcodes_stats']
    f_st = conf_f['fold_struct']

    target_list = [f"{inp_pref}{f_st}{bai_suff}", f"{inp_pref}{f_st}{gf_mat}", 
        f"{inp_pref}{f_st}{gf_features}", f"{inp_pref}{f_st}{barcodes_stats}"]
    return target_list


# Add support for multiome
# Add CollectInsertSizeMetrics for ATAC part of multiome
def targets_PICARD(conf_f, progs='all') -> list:
    
    inp_pref = conf_f['STARsolo_pipeline']['bams_dir']
    rna_seq_suff = conf_f['picard_pipeline']['rnaseq_metrics']
    gc_met_suff = conf_f['picard_pipeline']['gc_bias_metrics']
    gc_summ_suff = conf_f['picard_pipeline']['gc_summary_metrics']
    f_st=conf_f['fold_struct']

    files_dict = {
        'rnaseq':rna_seq_suff, 
        'gc':[gc_met_suff, gc_summ_suff],
        }
    target_list = []
    if progs == 'all':
       # Run both PICARD programs
        for k, val in files_dict.items():
            if isinstance(val, list):
                for l_val in val:
                    target_list.append(f"{inp_pref}{f_st}{l_val}")

            else:
                target_list.append(f"{inp_pref}{f_st}{val}")

    elif progs == "rnaseq":
        # Only RNAseq metrics
        target_list.extend( [f"{inp_pref}{f_st}{val}"  for val in files_dict[progs]] )

    elif progs == "gc":
        # Only GC bias metrics
        target_list.extend( [f"{inp_pref}{f_st}{val}"  for val in files_dict[progs]] )
    else:
        # print("Wrong input; Check for editing errors!")
        return []

    return target_list



# Temporary function, will change to targets_gt_demux
def targets_cellSNP(conf_f, progs=None, multiome=False) -> list:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    out_dir = conf_f['gt_demux_pipeline']['cellsnp_dir']
    fs_gt = conf_f['fold_struct_gt_demux']
    cells_suff = conf_f['gt_demux_pipeline']['cellsnp_cells']
    base_suff = conf_f['gt_demux_pipeline']['cellsnp_base']
    for d in sub_dir:
        # target_list = [f"{out_dir}{fs_gt}{cells_suff}", f"{out_dir}{fs_gt}{base_suff}"]
        target_list.extend([
            os.path.join(f"{out_dir}{fs_gt}", d, f"{cells_suff}"), 
            os.path.join(f"{out_dir}{fs_gt}", d, f"{base_suff}")
            ])

    # STARsolo* + PICARD (any) progs
    if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
        target_list.extend(targets_PICARD(conf_f=conf_f, progs=progs))

    return target_list



def targets_gt_demux(conf_f, progs=None, h5ad=False, multiome=False) -> list:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    if not h5ad:
        for d in sub_dir:
        
            out_dir = conf_f['gt_demux_pipeline']['vireosnp_dir']
            fs_gt = conf_f['fold_struct_gt_demux']
            suff = config['gt_demux_pipeline']['donors_classification']

            target_list.append(os.path.join(f"{out_dir}{fs_gt}", d ,f"{suff}"))
            # target_list = [f"{out_dir}{fs_gt}{suff}"]
    # For multiome combine cDNA and ATAC based demultiplexing into
    # a single h5ad
    else:
        out_dir = conf_f['gt_demux_pipeline']['final_count_matrix_dir']
        fs_gt = conf_f['fold_struct_demux']
        suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']

        # target_list = [f"{out_dir}{fs_gt}{suff}"]

        target_list.append(os.path.join(f"{out_dir}{fs_gt}{suff}"))


    # STARsolo* + PICARD (any) progs
    if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
        target_list.extend(targets_PICARD(conf_f=conf_f, progs=progs))

    return target_list


# For multi-vcf_inputs
def targets_gt_demux2(conf_f, progs=None, multiome=False) -> list:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    for d in sub_dir:
        out_dir = conf_f['gt_demux_pipeline']['final_count_matrix_dir']
        fs_gt = conf_f['fold_struct_demux']
        
        # target_list = [f"{out_dir}{fs_gt}" ]
        target_list.append(os.path.join(f"{out_dir}", d, f"{fs_gt}"))
    # STARsolo* + PICARD (any) progs
    if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
        target_list.extend(targets_PICARD(conf_f=conf_f, progs=progs))
    
    return target_list

# NEED TO WORK
def targets_SplitBams(conf_f, progs=None, multiome=False) -> list:
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    for d in sub_dir:
        if config['gt_check']:
            out_dir = conf_f['split_bams_pipeline']['split_bams_proxy_dir2']
        else:
            out_dir = conf_f['split_bams_pipeline']['split_bams_proxy_dir']
        op = conf_f['fold_struct_bam_split1']

        target_list.append(os.path.join(f"{out_dir}", f"{op}", d))
    
    # STARsolo* + PICARD (any) progs
    if progs == 'all' or progs == 'rnaseq' or progs == 'gc':
        target_list.extend(targets_PICARD(conf_f=conf_f, progs=progs, 
            ))
    
    return target_list


def targets_resolve_swaps_gt_demux(conf_f) -> list:
    out_dir = conf_f['finalize_demux']['out_dir']
    fs = conf_f['fold_struct_gt_demux_redo']
    dons_file = conf_f['finalize_demux']['donor_info_swap_file']
    suff = conf_f['gt_demux_pipeline']['donors_classification']
    
    return [f"{out_dir}{fs}{suff}"]


def targets_gt_demux_identify_swaps(conf_f, multiome=False) -> list:

    out_dir = conf_f['identify_swaps']['mbv_out_dir']
    fs = conf_f['fold_struct_swaps_check']
    suff = conf_f['identify_swaps']['mbv_suffix']
    sub_dir = ["ATAC", "cDNA"] if multiome else ['']
    target_list = []
    for d in sub_dir:    
        target_list.append(os.path.join(f"{out_dir}", d, f"{fs}{suff}"))

    return target_list


def targets_multibamsummary(conf_f) -> list:
    out_dir = conf_f['multiBamSummary']['outdir']
    fs = conf_f['fold_struct_deeptools']
    suff = ".npz"
    
    return [f"{out_dir}{fs}{suff}"]


def targets_multibamsummaryPlotCorr(conf_f) -> list:
    out_dir = conf_f['plotCorrelation']['outdir']
    fs = conf_f['fold_struct_deeptools']
    suff = conf_f['plotCorrelation']['out_fmt']
    
    return [f"{out_dir}{fs}{suff}"]


# Need to add support for PICARD metrics
# and create h5ad
def targets_multiome(conf_f, last, progs=None, 
    h5ad=False, multiome=False) -> list:
    out_dir = conf_f['cellranger_arc_count']['bams_dir']
    if last == 'alignment':
        return [
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/barcodes.tsv.gz",
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/features.tsv.gz",
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/matrix.mtx.gz",
        ]
    elif last == 'vireo':
        return targets_gt_demux(conf_f=conf_f, progs=None, 
            h5ad=h5ad, multiome=multiome)

    elif last == 'splitBams':
        return targets_SplitBams(conf_f=conf_f, progs=progs, multiome=multiome)

    elif last == 'identifySwaps':
        return targets_gt_demux_identify_swaps(conf_f=conf_f, multiome=multiome)


# To run STARsolo* + kb pipeline + (optional)PICARD progs
def targets_all(conf_f, progs=None) -> list:
    
    demuxed_mat_dir = conf_f['hashsolo_demux_pipeline']['final_count_matrix_dir']
    demuxed_info_dir = conf_f['hashsolo_demux_pipeline']['demultiplex_info_dir']    
    demux_mat_suff = conf_f['hashsolo_demux_pipeline']['final_count_matrix_h5ad']
    demux_info_suff = conf_f['hashsolo_demux_pipeline']['demultiplex_info']
    # f_st_bam = conf_f['fold_struct']
    f_st = conf_f['fold_struct_demux']

    # STARsolo* + PICARD (any) progs
    if progs is not None:
        target_list = targets_PICARD(conf_f=conf_f, progs=progs)

    else: # For any wrong values to progs or progs == None or just PICARD == False or everything else, run STARsolo*
        target_list = targets_STARsolo(conf_f=conf_f)


    # kb pipeline
    target_list.extend( [f"{demuxed_mat_dir}{f_st}{demux_mat_suff}", f"{demuxed_info_dir}{f_st}{demux_info_suff}"] )
    
    return target_list



# Add wildcards info here in the 'expand' function
def produce_targets(conf_f: pd.DataFrame, last_step: str, wc_d: dict) -> list:
    # global ONLY_SOLO, ONLY_VIREO, BOTH_DEMUX
    # conf_f['last_step'] is expected to be None if yaml file is provided for the config['select_fastqs']
    metrics = conf_f['picard_metrics'].lower() if conf_f['picard_metrics'] is not None else None
    
    
    if last_step is not None:
        target_step = last_step.lower()
        create_h5ad = conf_f['create_final_h5ad'] if 'gt_demux' in target_step else False
        multiome = True if 'multiome' in target_step else False
        # 'all' is only for the modules: STARsolo, PICARD (both) and calico_solo
        # and vireo (demux) *
        # In the future add create h5ad and split bams for gt check
        # NEEDS RE-CONSIDERATION
        if target_step == "all":
            target_files = targets_all(conf_f=conf_f, progs=metrics)
            target_files.extend(targets_gt_demux(conf_f=conf_f, progs=metrics))
            # global_vars.BOTH_DEMUX = True
            final_target_list = [expand(f"{target}", zip, **wc_d) for target in target_files] # if multiple wildcards are present then either 'zip' or 'product' (default)

        # NEEDS RE-CONSIDERATION
        elif target_step == "all_multi_vcf":
            target_files = targets_all(conf_f=conf_f, progs=metrics)
            # global_vars.BOTH_DEMUX = True
            final_target_list = [expand(f"{target}", zip, **wc_d) for target in target_files] # if multiple wildcards are present then either 'zip' or 'product' (default)
            # Adding gt_demux_multi_vcf_options
            target_files = targets_gt_demux2(conf_f=conf_f)
            suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']
            # global_vars.ONLY_VIREO = True # Set this or ADD_VIREO by figuring out or asking user input

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

                # temp_list= [expand(f"{target}", **wc_d) for target in target_files][0] # Single wildcard

            # Single vcf out of a set of multiple vcfs
            elif isinstance(VCF_TYPE, str):
                
                final_target_list.extend(
                    [expand(f"{target}_{VCF_TYPE}{suff}", zip, **wc_d) for target in target_files]
                    ) # Multiple wildcards example
                # final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", **wc_d) for target in target_files] # Single wildcard
            else:
                raise ValueError("The wildcard used to test multi_vcf "
                    "input (VCF_TYPE) is of unexpected type! Please check")

        elif target_step == "starsolo_kb_solo":
            # global_vars.ONLY_SOLO = True
            target_files = targets_all(conf_f=conf_f, progs=metrics)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        elif target_step == "starsolo":
            if metrics is None:
                target_files = targets_STARsolo(conf_f=conf_f)
            else:
                target_files = targets_PICARD(conf_f=conf_f, progs=metrics)

            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]
            
        # elif target_step == "starsolo_picard":
        #     target_files = targets_PICARD(conf_f=conf_f, progs='all')

        #     final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        elif target_step == "starsolo_cellsnp":
            target_files = targets_cellSNP(conf_f=conf_f, progs=metrics)

            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]    

            
        elif target_step == "starsolo_gt_demux":
            target_files = targets_gt_demux(conf_f=conf_f, progs=metrics, 
                            h5ad=create_h5ad, multiome=multiome)
            # global_vars.ONLY_VIREO = True
            # if not create_h5ad:
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]
            # else:
            #     final_target_list = target_files

        # Multi_vcf inputs or when a subdir with the vcf name is needed
        # NEEDS RE-CONSIDERATION
        elif target_step == "starsolo_gt_demux_multi_vcf":
            target_files = targets_gt_demux2(conf_f=conf_f)
            suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']
            # global_vars.ONLY_VIREO = True # Set this or ADD_VIREO by figuring out or asking user input

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

        elif target_step == "starsolo_split_bams":
            target_files = targets_SplitBams(conf_f=conf_f, progs=metrics)

            # global_vars.ONLY_SOLO = True
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        elif target_step == "starsolo_split_bams_gt_demux":
            target_files = targets_SplitBams(conf_f=conf_f, progs=metrics)
            suff = ".txt"
            # global_vars.ONLY_VIREO = True
            # temp_list= [expand(f"{target}", zip, num=round_num, **wc_d) for target in target_files][0] # Multiple wildcards example
            # final_target_list= [expand(f"{target}{suff}", zip, **wc_d) for target in target_files][0] # Single wildcard
            
            final_target_list = []
            for id, target in enumerate(target_files):
                if id == 0:
                    final_target_list.extend(expand(f"{target}{suff}", zip, **wc_d))
                else:
                    final_target_list.extend(expand(f"{target}", zip, **wc_d))

            # final_target_list.append(f"{config['split_bams_pipeline_gt_demux']['ngscheckmate_res_dir']}/{config['split_bams_pipeline_gt_demux']['proj_name']}_all.txt")

        # Multi_vcf inputs or when a subdir with the vcf name is needed
        elif target_step == "starsolo_split_bams_gt_demux_multi_vcf":
            # Add final steps for split_bams
            target_files = targets_SplitBams(conf_f=conf_f, progs=metrics)
            suff = ".txt"
            # If multiple vcfs per sample needs to be run or not
            if isinstance(VCF_TYPE, list):
                
                temp_list = [expand(f"{target}", zip, **wc_d) for target in target_files][0] # Multiple wildcards example
                # temp_list= [expand(f"{target}", **wc_d) for target in target_files][0] # Single wildcard
                final_target_list = [expand(f"{t}_{{vcf_type}}{suff}", vcf_type=VCF_TYPE) for t in temp_list]

            # Single vcf out of a set of multiple vcfs
            elif isinstance(VCF_TYPE, str):
                
                final_target_list = [expand(f"{target}_{VCF_TYPE}{suff}", zip, **wc_d) for target in target_files] # Multiple wildcards example
                # final_target_list= [expand(f"{target}_{VCF_TYPE}{suff}", **wc_d) for target in target_files] # Single wildcard
            else:
                raise ValueError("The wildcard used to test multi_vcf input (VCF_TYPE) is of unexpected type! Please check")

            # Add final steps for gt_demux
            # target_files = targets_gt_demux2(conf_f=conf_f, progs=metrics)
            # suff = config['gt_demux_pipeline']['final_count_matrix_h5ad']
            # ONLY_VIREO = True # Set this or ADD_VIREO by figuring out or asking user input

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

            # target_files = targets_PICARD(conf_f=conf_f, progs=metrics)
            # final_target_list.extend([expand(f"{target}", zip, **wc_d) for target in target_files])
            # final_target_list.extend(final_target_list2)
            # final_target_list.append(f"{config['split_bams_pipeline_gt_demux']['ngscheckmate_res_dir']}/{config['split_bams_pipeline_gt_demux']['proj_name']}_all.txt")

        elif target_step == "starsolo_gt_demux_identify_swaps":
            target_files = targets_gt_demux_identify_swaps(conf_f=conf_f)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        elif target_step == "starsolo_resolve_swaps_gt_demux":
            target_files = targets_resolve_swaps_gt_demux(conf_f=conf_f)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        elif target_step == "multibamsummaryplotcorr":
            # global_vars.ONLY_SOLO = True
            target_files = targets_multibamsummaryPlotCorr(conf_f=conf_f)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files[0]]

        elif target_step == "multibamsummary":
            # global_vars.ONLY_SOLO = True
            target_files = targets_multibamsummary(conf_f=conf_f)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files[0]]

        elif target_step == "multiome_alignment":
            target_files = targets_multiome(conf_f=conf_f, last="alignment")
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        # Support creation of h5ad or not
        elif "multiome_gt_demux" in target_step:
            target_files = targets_multiome(conf_f=conf_f, last="vireo", progs=metrics, 
            multiome=multiome, h5ad=create_h5ad)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        elif "multiome_split_bams_gt_demux" in target_step:
            target_files = targets_multiome(conf_f=conf_f, last="splitBams", progs=metrics,
            multiome=multiome, h5ad=create_h5ad)
            suff = ".txt"
            # global_vars.ONLY_VIREO = True
            # temp_list= [expand(f"{target}", zip, num=round_num, **wc_d) for target in target_files][0] # Multiple wildcards example
            # final_target_list= [expand(f"{target}{suff}", zip, **wc_d) for target in target_files][0] # Single wildcard
            
            final_target_list = []

            for id, target in enumerate(target_files):
                if id <= 1:
                    final_target_list.extend(expand(f"{target}{suff}", zip, **wc_d))
                else:
                    final_target_list.extend(expand(f"{target}", zip, **wc_d))

        elif target_step == "multiome_gt_demux_identify_swaps":
            target_files = targets_multiome(conf_f=conf_f, last="identifySwaps", progs=metrics,
            multiome=multiome, h5ad=create_h5ad)
            final_target_list= [expand(f"{target}", zip, **wc_d) for target in target_files]

        else:
            # print("Wrong inputs to produce_targets function!!")
            return []


        # Fail-safe return statement
        return final_target_list


    else:
        raise ValueError("The variable last_step is set to None!")



