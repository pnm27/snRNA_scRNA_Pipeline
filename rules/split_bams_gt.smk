# For splitting bams
# def get_donors(wildcards, input):
# 	hash_file=pd.read_csv(input.barcodes_vs_donor, sep='\t')
# 	# Make this generalized with loc - bc hash co0lumn names in config.yaml file
# 	return ' '.join(hash_file.iloc[:, 0].unique().tolist())


def subset_to_chr(wildcards):
    if config['chr_prefix'] == None or config['split_bams_pipeline_gt_demux']['subset_chr'].startswith(config['chr_prefix']):
        return config['split_bams_pipeline_gt_demux']['subset_chr']
    else:
        return config['chr_prefix'] + config['split_bams_pipeline_gt_demux']['subset_chr'] 


def get_hash_file(wildcards):
    # Not updated for multi_modules
    if config['last_step'] is not None and config['last_step'].lower() == "starsolo_split_bams":
        return f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt" # with headers
    # Create hash from vireo_SNPs output
    # If donor names aren't the same then use a converter file to do so
    # The output file should have headers (the same as when calico_solo's output when run through the rule create_inp_splitBams)
    elif config['last_step'] is not None and config['last_step'].lower() == "starsolo_split_bams_gt_demux": # change this soon
        return f"{config['split_bams_pipeline_gt_demux']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"
    else:
        raise ValueError("Something is wrong with the expected input for the rule split_bams")


def get_bam_to_split(wildcards):
    if config['split_bams_pipeline_gt_demux']['subset_chr'] is None:
        return f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"
    else:
        return f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['short_bam']}"


def get_mito(wildcards):
    if config['chr_prefix'] == None or config['mito'].startswith(config['chr_prefix']):
        return config['mito']
    else:
        return config['chr_prefix'] + config['mito'] 



rule create_inp_splitBams:
    input:
        f"{config['demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['demux_pipeline']['final_count_matrix_h5ad']}"

    output:
        f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"

    params:
        overwrite=config['split_bams_pipeline']['overwrite']


    threads: 1

    resources:
        mem_mb=allocate_mem_CIS,
        time_min=10

    conda: "../envs/basic_sctools.yaml"

    shell:
        """
        if [ "{params.overwrite}" == True ]; then
            python3 helper_py_scripts/create_inp_splitBam.py {input} {output} --overwrite

        else
            python3 helper_py_scripts/create_inp_splitBam.py {input} {output}

        fi
        sleep 60
        """



rule create_inp_splitBams_gt_demux:
    input:
        f"{config['gt_demux_pipeline']['vireosnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['donors_classification']}"

    output:
        f"{config['split_bams_pipeline_gt_demux']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"

    params:
        conv_df=config['gt_conv']

    threads: 1

    resources:
        mem_mb=allocate_mem_CISPD,
        time_min=10

    conda: "../envs/basic_sctools.yaml"

    shell: 
        """
        if [[ "{params.conv_df}" != "None" ]]; then 
            python3 helper_py_scripts/create_inp_splitBam_gt_demux.py \
            {input[0]} {output[0]} --converter_file {params.conv_df}
        else
            python3 helper_py_scripts/create_inp_splitBam_gt_demux.py \
            {input[0]} {output[0]}
        fi
        """




rule filt_chr_bams:
    input:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['short_bam']}" # generalize this

    params:
        temp_bam=f"{config['split_bams_pipeline']['sort_temp_dir']}temp_bam/{{id1}}.bam",
        sub_chr=subset_to_chr


    threads: 1

    resources:
        mem_mb=allocate_mem_FCB,
        time_min=allocate_time_FCB

    conda: "../envs/pysam.yaml"

    envmodules:
        "samtools"
        
    shell:
        """
        ml samtools
        set -x
        samtools view {input} {params.sub_chr} -bho {output}
        set +x
        """


rule split_bams:
    input:
        short_bam=get_bam_to_split,
        barcodes_vs_donor=get_hash_file

    params:
        # split_at=config['split_bams_pipeline']['bc_per_donor'], # Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
        # temp_dir=f"{config['split_bams_pipeline']['temp_dir']}",
        temp_bam_per_cell_dir=f"{config['split_bams_pipeline_gt_demux']['new_temp_dir2']}{config['fold_struct_bam_split2']}",
        split_bams_dir=f"{config['split_bams_pipeline_gt_demux']['split_bams_dir2']}{config['fold_struct_bam_split2']}",
        per_donor_log_dir=config['split_bams_pipeline_gt_demux']['per_donor_split_log_dir2'],
        time_limit_per_donor=config['split_bams_pipeline_gt_demux']['time_per_donor']
        # n_donors=get_donors


    threads: 1

    resources:
        mem_mb=allocate_mem_SB,
        time_min=10

    output:
        f"{config['split_bams_pipeline_gt_demux']['split_bams_proxy_dir2']}{config['fold_struct_bam_split1']}.txt"  # Proxy to the output

    run:
        hash_file=pd.read_csv(input.barcodes_vs_donor, sep='\t')
        # bam_files=0
        proc_donors=[d for d in hash_file.iloc[:, 0].unique() if not os.path.isfile(os.path.join(params.split_bams_dir, f"{d}.bam"))]
        job_name_l=[]
        samp_num=('-'.join(os.path.basename(input.short_bam).replace('_chr1.bam', '').split('-')[1:3])) # generalize the name
        if not os.path.isdir(params.temp_bam_per_cell_dir):
            os.makedirs(params.temp_bam_per_cell_dir)

        if not os.path.isdir(params.split_bams_dir):
            os.makedirs(params.split_bams_dir)

        if not os.path.isdir(params.per_donor_log_dir):
            os.makedirs(params.per_donor_log_dir)

        # If proc_donors is empty don't run
        if proc_donors:
            for donor in proc_donors:
                # bam_files+=1
                # Change this nameing style
                jname=f"NPSAD_{samp_num}_" + list(shell("uuidgen", iterable=True))[0]
                job_name_l.append(jname)
                shell("""
            		jid=$(bsub -J {j} -P acc_CommonMind -q express -n 1 -R span[hosts=1] -R rusage[mem=200] -W {params.time_limit_per_donor} -oo {params.per_donor_log_dir}{j}.stdout -eo {params.per_donor_log_dir}{j}.stderr -L /bin/bash "bash helper_sh_scripts/create_per_donor_bams_gt.bash {i} {input.barcodes_vs_donor} {params.temp_bam_per_cell_dir} {params.split_bams_dir} {input.short_bam}")
                    jid=$(echo $jid | head -n1 | cut -d '<' -f2 | cut -d '>' -f1)
                    echo "Submitted script for donor {i} with jobid: ${{jid}}" >> {output}
                    sleep 10
                """, i=donor, j=jname)

            with open(output[0], "a") as fout:
                fout.write(f"Number of 'new' bam files expected at the completion of all the scripts for the bam file {input.short_bam} is {len(job_name_l)}\n")


        else:
            with open(output[0], "a") as fout:
                fout.write(f"Skipped bam file {input.short_bam} as all {len(hash_file.iloc[:, 0].unique())} donor file(s) was(were) already present in the given output_folder {params.split_bams_dir}")


# Integrate this later by adding checkpoints and also create
# The input file required
# Make sure this runs after each type of demux i.e. after calico_solo and after vireo
# NGScheckmate git has issues with indentations in ncm.py
# Open it up in an editor to correct them
# rule run_NGSCheckmate:
#     params:
#         # split_bams_dir=f"{config['split_bams_pipeline_gt_demux']['split_bams_dir2']}{config['fold_struct_bam_split2']}",
#         inp_file_list=config['split_bams_pipeline_gt_demux']['inp_file'],
#         ngscheckmate_dir=config['split_bams_pipeline_gt_demux']['ngscheckmate_dir'],
#         ref_fasta=config['genome_fasta'],
#         ref_bed=config['split_bams_pipeline_gt_demux']['ngscheckmate_ref_bed'],
#         inp_type=config['split_bams_pipeline_gt_demux']['ngs_inp_type'],
#         proj_name=config['split_bams_pipeline_gt_demux']['proj_name'],
#         out_dir={config['split_bams_pipeline_gt_demux']['ngscheckmate_res_dir']}
#         # n_donors=get_donors

#     threads: 5

#     resources:
#         mem_mb=allocate_mem_RN
#         # time_min=allocate_time_RN

#     output:
#         f"{config['split_bams_pipeline_gt_demux']['ngscheckmate_res_dir']}/{config['split_bams_pipeline_gt_demux']['proj_name']}_all.txt"
#     shell:
#         """
#         ml samtools
#         ml bcftools/1.15.1
#         export NCM_HOME={params.ngscheckmate_dir}NGSCheckMate
#         echo "REF={params.ref_fasta}" > ${{NCM_HOME}}/ncm.conf
#         echo "SAMTOOLS=samtools" >> ${{NCM_HOME}}/ncm.conf
#         echo "BCFTOOLS=bcftools" >> ${{NCM_HOME}}/ncm.conf
#         python ${{NCM_HOME}}/ncm.py -{params.inp_type} -l {params.inp_file_list} -bed {params.ref_bed} -N {params.proj_name} -O {params.out_dir}
#         sleep 30
#         """


# rule consolidate_gt_check_results:
#     input:
#         short_bam=get_bam_to_split,
#         barcodes_vs_donor=get_hash_file

#     params:
#         # split_at=config['split_bams_pipeline']['bc_per_donor'], # Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
#         # temp_dir=f"{config['split_bams_pipeline']['temp_dir']}",
#         temp_bam_per_cell_dir=f"{config['split_bams_pipeline_gt_demux']['new_temp_dir2']}{config['fold_struct_bam_split2']}",
#         split_bams_dir=f"{config['split_bams_pipeline_gt_demux']['split_bams_dir2']}{config['fold_struct_bam_split2']}",
#         per_donor_log_dir=config['split_bams_pipeline_gt_demux']['per_donor_split_log_dir2'],
#         time_limit_per_donor=config['split_bams_pipeline_gt_demux']['time_per_donor']
#         # n_donors=get_donors


#     threads: 1

#     resources:
#         mem_mb=allocate_mem_SB,
#         time_min=10

#     output:
#         f"{config['split_bams_pipeline_gt_demux']['split_bams_proxy_dir2']}{config['fold_struct_bam_split1']}.txt"  # Proxy to the output

#     run: