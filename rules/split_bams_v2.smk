import os, pandas as pd

def get_chr_pref(wildcards):
    if config['chr_prefix'] != None:
        return ''
    else:
        return config['chr_prefix']


def get_mito(wildcards):
    if config['chr_prefix'] == None or config['mito'].startswith(config['chr_prefix']):
        ret_val = config['mito'] if not config['mito'].endswith('-') else config['mito'][:-1]
        return ret_val
    else:
        ret_val = config['chr_prefix'] + config['mito']
        ret_val =  ret_val if not ret_val.endswith('-') else ret_val[:-1]
        return ret_val


# defining local rules
# localrules: create_bed

def get_donor_bams(wildcards):
    checkpoint_output = checkpoints.split_bams.get(**wildcards).output[0]
    temp_df = pd.read_csv(config['demultiplex']['final_demultiplex_info'], sep='\t')
    sample_name = wildcards.id1 # CHANGE THIS TO THE WILDCARD THAT CONTAINS THE SAMPLE NAME
    pool_col = config['demultiplex']['final_demultiplex_info_poolCol']
    donor_col = config['demultiplex']['final_demultiplex_info_donorCol']
    don_sep = config['demultiplex']['final_demultiplex_info_donorSep']
    # When all donors present in a pool are present in one row as:
    # pool_col  donor_col
    # poolA donor1, donor2, donor3, donor4
    # poolB donor5, donor6, donor7, donor8
    # ..... ......
    # etc.
    if don_sep is not None:
        donors = temp_df.loc[temp_df[pool_col] == sample_name, donor_col].values[0]
        return list(map(lambda x: checkpoint_output + x + ".bam", donors.split(don_sep)))
    else:
        donors = temp_df.loc[temp_df[pool_col] == sample_name, donor_col].values.tolist()
        return list(map(lambda x: checkpoint_output + x + ".bam", don_sep))


rule create_inp_splitBams:
    input:
        f"{config['demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['demux_pipeline']['final_count_matrix_h5ad']}"

    priority: 7

    output:
        f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"

    params:
        overwrite=config['split_bams_pipeline']['overwrite']


    threads: 2

    resources:
        mem_mb=1000,
        time_min=10

    shell:
        """
        if [ "{params.overwrite}" == True ]; then
            python3 helper_py_scripts/create_inp_splitBam.py {input} {output} --overwrite --split_by SubID_vs_w_gt --demux_suffix vs

        else
            python3 helper_py_scripts/create_inp_splitBam.py {input} {output}

        fi
        sleep 60
        """


rule bamfilt_by_CB:
    input:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}",
        f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"

    output:
        temp(f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}")
        

    params:
        temp_bc=f"{config['split_bams_pipeline']['sort_temp_dir']}{{num}}_{{id1}}_bc.txt"


    threads: 1

    resources:
        mem_mb=250,
        time_min=300

    shell:
        """
        mkdir -p {config[split_bams_pipeline][sort_temp_dir]}
        ml samtools
        cut -f2 <(tail -n +2 {input[1]}) > {params.temp_bc}
        samtools view -q 255 -D CB:{params.temp_bc} {input[0]} -bho {output}
        sleep 60
        samtools index {output}
        rm {params.temp_bc}
        sleep 10
        """


rule create_bed:
    input:
        config['genome_fasta']

    params:
        chr_prefix=get_chr_pref

    output:
        config['reg_chr_bed']

    shell:
        """
        if [ ! -f "{output}" ] ; then 
            awk -v pref="^{params.chr_prefix}[0-9XY]*\\t" "" $0 ~ pref {{printf("%s\\t0\\t%s\\n",\$1,\$2);}}"" {input} > {output}
        fi
        sleep 10
        """


# rule collect_mito_info:
#     input:
#         config['reg_chr_bed'],
#         filt_bam=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}"

#     params:
#         # split_at=config['split_bams_pipeline']['bc_per_donor'], # Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
#         # temp_dir=f"{config['split_bams_pipeline']['temp_dir']}",
#         temp_bam=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}{{num}}_{{id1}}.bam",
#         chr_mito=get_mito

#     output:
#         f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['mito_reads_file']}"

#     shell:
#         """
#             ml samtools
#             op={output}
#             if [ ! -d "${{op%$(basename ${{op}})}}" ]; then mkdir -p ${{op%$(basename ${{op}})}}; fi
#             if [ ! -f "{output}" ] || ! grep -E -q "Number of Mitochondrial reads after retaining only mapped reads for pooled sample is: [0-9]+" "${output}"
#             then
#                 echo -n "Number of Mitochondrial reads after retaining only mapped reads for pooled sample is: " > {output} && samtools view -c {input.filt_bam} {params.chr_mito} >> {output}
#             fi
#             if [ ! -f "{params.temp_bam}" ]; then samtools view -L {input[0]} -o {params.temp_bam} {input.filt_bam}; fi
#             sleep 20
#         """




checkpoint split_bams:
    input:
        config['reg_chr_bed'],
        filt_bam=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}",
        barcodes_vs_donor=f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt" # with headers

    priority: 7

    params:
        # split_at=config['split_bams_pipeline']['bc_per_donor'], # Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
        # temp_dir=f"{config['split_bams_pipeline']['temp_dir']}",
        temp_bam=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}{{num}}_{{id1}}.bam",
        temp_bam_per_cell_dir=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}",
        split_bams_dir=f"{config['split_bams_pipeline']['split_bams_dir']}{config['fold_struct_bam_split2']}",
        per_donor_log_dir=config['split_bams_pipeline']['per_donor_split_log_dir'],
        time_limit_per_donor=config['split_bams_pipeline']['time_per_donor'],
        chr_mito=get_mito

    threads: 1

    resources:
        mem_mb=allocate_mem_SB,
        time_min=30

    output:
        # f"{config['split_bams_pipeline']['split_bams_proxy_dir']}{config['fold_struct_bam_split1']}.txt",  # Proxy to the output
        # f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['mito_reads_file']}"
        directory(f"{config['split_bams_pipeline']['split_bams_dir']}{config['fold_struct_bam_split2']}")

    run:
        hash_file=pd.read_csv(input.barcodes_vs_donor, sep='\t')
        proc_donors=[d for d in hash_file.iloc[:, 0].unique() if not os.path.isfile(os.path.join(params.split_bams_dir, f"{d}.bam.bai"))]
        job_name_l=[]
        samp_num=('-'.join(os.path.basename(input.filt_bam).replace('_filt.bam', '').split('-')[1:3]))
        
        if not os.path.isdir(params.temp_bam_per_cell_dir):
            os.makedirs(params.temp_bam_per_cell_dir)

        if not os.path.isdir(params.split_bams_dir):
            os.makedirs(params.split_bams_dir)

        if not os.path.isdir(params.per_donor_log_dir):
            os.makedirs(params.per_donor_log_dir)

        # If proc_donors is not empty then run
        if proc_donors:
            with open(output[0], "a") as fout:
                fout.write("Going to process donor(s): {}".format(','.join(proc_donors)))

            for donor in proc_donors:
                jname=f"NPSAD_{samp_num}_" + list(shell("uuidgen", iterable=True))[0]
                job_name_l.append(jname)
                shell("""
                    if [ ! -d "{params.per_donor_log_dir}" ]; then mkdir -p {params.per_donor_log_dir}; fi
                    jid=$(bsub -J {j} -P acc_CommonMind -q express -n 1 -R span[hosts=1] -R rusage[mem=200] -W {params.time_limit_per_donor} -oo {params.per_donor_log_dir}{j}.stdout -eo {params.per_donor_log_dir}{j}.stderr -L /bin/bash "bash helper_sh_scripts/create_per_donor_bams.bash {i} {input.barcodes_vs_donor} {params.temp_bam_per_cell_dir} {params.split_bams_dir} {input.filt_bam} {output[1]} {input[0]} {params.chr_mito}")
                    jid=$(echo $jid | head -n1 | cut -d '<' -f2 | cut -d '>' -f1)
                    echo "Submitted script for donor {i} with jobid: ${{jid}}" >> {output[0]}
                    sleep 10
                """, i=donor, j=jname)

            with open(output[0], "a") as fout:
                fout.write(f"Number of 'new' bam files expected at the completion of all the scripts for the bam file {input.filt_bam} is {len(job_name_l)}")

        else:
            with open(output[0], "a") as fout:
                fout.write(f"Skipped bam file {input.filt_bam} as all {len(hash_file.iloc[:, 0].unique())} donor file(s) was(were) already present in the given output_folder {params.split_bams_dir}")


rule delete_filt_bam:
    input:
        get_donor_bams

    params:
        filt_bam=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}"
    output:
        f"{config['split_bams_pipeline']['split_bams_proxy_dir']}{config['fold_struct_bam_split1']}.txt"

    shell:
        """
            op={output}
            if [ ! -d "${{op%$(basename ${{op}})}}" ]; then mkdir -p ${{op%$(basename ${{op}})}}; fi
            echo 'Created all individual bam files: ' >>> {output}
            echo '{input}' >>> {output}
            rm {params.filt_bam} && echo 'Deleted temporary filtered file {{params.filt_bam}}' >>> {output} || echo 'Could not delete temporary file {{params.filt_bam}}' >>> {output}
        """