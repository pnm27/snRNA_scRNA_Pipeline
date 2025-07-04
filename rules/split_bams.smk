import os, pandas as pd


def get_inp_splitBam(wildcards):
    if config['split_bams_pipeline']['split_by']['input'].lower() == 'raw':
        if config['split_bams_pipeline']['split_by']['demux'].lower() in ['vireo', 'vs']:
            return f"{config['gt_demux_pipeline']['vireosnp_dir']}{config['fold_struct_gt_demux']}{config['gt_demux_pipeline']['donors_classification']}"
        elif config['split_bams_pipeline']['split_by']['demux'].lower() in \
            ["cs", "calico", "calico_solo", "hashsolo"]:
            return f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
    else:
        # if len(config['split_bams_pipeline']['split_by']['demux']) == 1:
        #     if config['split_bams_pipeline']['split_by']['demux'][0].lower() in ['vireo', 'vs'] or \
        #       config['split_bams_pipeline']['split_by']['demux'][0].lower() in \
        #       ["cs", "calico", "calico_solo", "hashsolo"]:
        #         return f"{config['demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['demux_pipeline']['final_count_matrix_h5ad']}"
        # else:
        return f"{config['hashsolo_demux_pipeline']['final_count_matrix_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['final_count_matrix_h5ad']}"

def get_bam(wildcards):
    if 'multiome' in config['last_step'].lower():
        if 'cdna' in wildcards.pool.lower(): # WILDCARDS
            return (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
                
        else:
            return (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
    else:
        return (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"
        )


def get_h5ad_cols(wildcards):
    return ' '.join()


def get_chr_pref(wildcards):
    if config['chr_prefix'] is None:
        return ''
    else:
        return config['chr_prefix']


def get_mito(wildcards):
    if config['chr_prefix'] == None or config['mito_prefix'].startswith(config['chr_prefix']):
        ret_val = config['mito_prefix'] if not config['mito_prefix'].endswith('-') else config['mito_prefix'][:-1]
        return ret_val
    else:
        ret_val = config['chr_prefix'] + config['mito_prefix']
        ret_val =  ret_val if not ret_val.endswith('-') else ret_val[:-1]
        return ret_val


def subset_to_chr(wildcards):
    if config['chr_prefix'] == None or str(config['split_bams_pipeline']['subset_chr']).startswith(config['chr_prefix']):
        return config['split_bams_pipeline']['subset_chr']
    else:
        return config['chr_prefix'] + str(config['split_bams_pipeline']['subset_chr'])


def get_bam_to_split(wildcards):
    fullBam = ""
    subBam = ""
    filtBam = ""
    if 'multiome' in config['last_step'].lower():
        if 'cdna' in wildcards.pool.lower(): # WILDCARDS
            fullBam += (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
            ).format(pool=wildcards.pool.split('/')[0])
            subBam += (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]).replace('.bam', config['split_bams_pipeline']['short_bam']) # WILDCARDS
            filtBam += (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['gex_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]).replace('.bam', config['split_bams_pipeline']['filt_bam']) # WILDCARDS
        else:
            fullBam += (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]) # WILDCARDS
            subBam += (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]).replace('.bam', config['split_bams_pipeline']['short_bam']) # WILDCARDS
            filtBam += (
                f"{config['cellranger_arc_count']['bams_dir']}"
                f"{{pool}}/{config['cellranger_arc_count']['atac_bam']}"
            ).format(pool=wildcards.pool.split('/')[0]).replace('.bam', config['split_bams_pipeline']['filt_bam']) # WILDCARDS
    else:
        fullBam += (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"
            )
        subBam += (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['split_bams_pipeline']['short_bam']}"
        )
        filtBam += (
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}"
        )
        
    if config['gt_check']:
        if config['split_bams_pipeline']['subset_chr'] is None:
            return fullBam
        else:
            return subBam
    else:
        return filtBam


def get_mito_file(wildcards):
    if config['gt_check']:
        return False
    elif not config['gt_check']:
        return f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['mito_reads_file']}"


# Resource Allocation ------------------
def allocate_mem_FCB(wildcards, attempt):
    return 150*attempt+150


def allocate_time_FCB(wildcards, attempt):
    return 20*attempt+20


def allocate_mem_SB(wildcards, attempt):
    return 75*attempt+75

# --------------------------------------

# FIX THE INPUT: INPUT SHOULD TAKE IN EITHER VIREO OR CALICO_SOLO'S H5AD AS INPUT, AS REQUIRED.
rule create_inp_splitBams:
    input:
        get_inp_splitBam

    # priority: 7

    output:
        f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"

    params:
        # overwrite=config['split_bams_pipeline']['overwrite']
        conv=config['split_bams_pipeline']['donor_name_converter']['file'],
        from_col=config['split_bams_pipeline']['donor_name_converter']['from_column'],
        to_col=config['split_bams_pipeline']['donor_name_converter']['to_column'],
        demux_method=config['split_bams_pipeline']['split_by']['demux'],
        inp_ext=config['split_bams_pipeline']['split_by']['input'],
        # n_methods=len(config['split_bams_pipeline']['split_by']['demux']),
        # demux_suffixes=len(config['split_bams_pipeline']['split_by']['suffix']),
        h5ad_col=config['split_bams_pipeline']['split_by']['column']

    # For snakemake < v8
    # threads: 2

    resources:
        cpus_per_task=2, # For snakemake > v8
        mem_mb=1000,
        time_min=10

    shell:
        """
        cmd_str="{input} {output} --split_by {params.demux_method} "
        if [[ "{params.conv}" != "None" && "{params.conv}" != "False" ]]; then
            cmd_str+="--conv_file {params.conv} "
            cmd_str+="--conv_file_from_col {params.from_col} "
            cmd_str+="--conv_file_to_col {params.to_col} "
            if [[ "{params.h5ad_col}" != "None" && "{params.inp_ext}" == "h5ad" ]]; then
                cmd_str+="--h5ad_donor_column {params.h5ad_col} "
            fi
        fi

        python3 helper_py_scripts/create_inp_splitBam_consolidated.py ${{cmd_str}}

        sleep 60
        """

# Generalize Output
rule bamfilt_by_CB:
    input:
        bam=get_bam,
        hash_file=f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt"

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['filt_bam']}"
        

    params:
        # NUM = pool, ID1=donor
        temp_bc=lambda wc: f"{config['split_bams_pipeline']['sort_temp_dir']}{{pool}}_bc.txt"


    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=250,
        time_min=300

    conda: "../envs/pysam.yaml"

    envmodules:
        "samtools"

    shell:
        """
        mkdir -p {config[split_bams_pipeline][sort_temp_dir]}
        cut -f2 <(tail -n +2 {input.hash_file}) > {params.temp_bc}
        samtools view -q 255 -D CB:{params.temp_bc} {input.bam} -bho {output}
        samtools index {output}
        rm {params.temp_bc}
        sleep 10
        """


rule filt_chr_bams:
    input:
        get_bam

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['split_bams_pipeline']['short_bam']}" # generalize this

    params:
        sub_chr=subset_to_chr

    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_FCB,
        time_min=allocate_time_FCB

    conda: "../envs/pysam.yaml"

    envmodules:
        "samtools"
        
    shell:
        """
        samtools view {input} {params.sub_chr} -bho {output}
        samtools index {output}
        """


rule filt_chr_bams_multiome:
    input:
        get_bam

    output:
        f"{config['cellranger_arc_count']['bams_dir']}{config['fold_struct']}{{bam}}{config['split_bams_pipeline']['short_bam']}" # generalize this

    params:
        sub_chr=subset_to_chr

    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_FCB,
        time_min=allocate_time_FCB

    conda: "../envs/pysam.yaml"

    envmodules:
        "samtools"

    shell:
        """
        samtools view {input} {params.sub_chr} -bho {output}
        samtools index {output}
        """


rule create_bed:
    input:
        config['genome_fasta']

    params:
        chr_prefix=get_chr_pref

    output:
        config['reg_chr_bed']

    shell:
        r"""
        if [ ! -f "{output}" ] ; then 
            grep -E "^>{params.chr_prefix}[0-9]+|X|Y|MT" "{input}" | \
            awk 'BEGIN{{OFS="\t"}}{{split($3,a,":");gsub(">", "", $1);printf("%s\\t0\\t%s\\n",\$1,a[5]);}}' > "{output}"
        fi
        sleep 10
        """


rule split_bams:
    input:
        filt_bam=get_bam_to_split,
        barcodes_vs_donor=f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt" # with headers

    params:
        # split_at=config['split_bams_pipeline']['bc_per_donor'], # Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
        # temp_dir=f"{config['split_bams_pipeline']['temp_dir']}",
        # temp_bam=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}{{pool}}_{{donor}}.bam", # wildcards
        temp_bam_per_cell_dir=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}",
        split_bams_dir=f"{config['split_bams_pipeline']['split_bams_dir']}{config['fold_struct_bam_split2']}",
        per_donor_log_dir=config['split_bams_pipeline']['per_donor_split_log_dir'],
        time_limit_per_donor=config['split_bams_pipeline']['time_per_donor'],
        chr_mito=get_mito,
        gt_check=config['gt_check'],
        mito_info=get_mito_file,
        filter_bed=None if config['gt_check'] else config['reg_chr_bed']

    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_SB,
        time_min=30

    output:
        f"{config['split_bams_pipeline']['split_bams_proxy_dir']}{config['fold_struct_bam_split1']}.txt"  # Proxy to the output

    run:
        hash_file=pd.read_csv(input.barcodes_vs_donor, sep='\t')
        proc_donors=[d for d in hash_file.iloc[:, 0].unique() if not os.path.isfile(os.path.join(params.split_bams_dir, f"{d}.bam.bai"))]
        job_name_l=[]
        samp_num=('-'.join(os.path.basename(input.filt_bam).replace('_filt.bam', '').split('-')[1:3])) # project-specific
        # If proc_donors is not empty then run
        if proc_donors:
            with open(output[0], "a") as fout:
                fout.write("Going to process donor(s): {}\n".format(','.join(proc_donors)))

            if params.gt_check:
                for donor in proc_donors:
                    jname=f"{samp_num}_" + list(shell("uuidgen", iterable=True))[0]
                    job_name_l.append(jname)
                    shell("""
                        if [ ! -d "{params.per_donor_log_dir}" ]; then mkdir -p {params.per_donor_log_dir}; fi
                        jid=$(bsub -J {j} -P acc_CommonMind -q express -n 1 -R span[hosts=1] \
                            -R rusage[mem=200] -W {params.time_limit_per_donor} \
                            -oo {params.per_donor_log_dir}{j}.stdout -eo {params.per_donor_log_dir}{j}.stderr \
                            -L /bin/bash "bash helper_sh_scripts/create_per_donor_bams.bash {i} \
                            {input.barcodes_vs_donor} {params.temp_bam_per_cell_dir} \
                            {params.split_bams_dir} {input.filt_bam} ")
                        jid=$(echo $jid | head -n1 | cut -d '<' -f2 | cut -d '>' -f1)
                        echo "Submitted script for donor {i} with jobid: ${{jid}}" >> {output[0]}
                        sleep 10
                    """, i=donor, j=jname)
            else:
                for donor in proc_donors:
                    jname=f"{samp_num}_" + list(shell("uuidgen", iterable=True))[0]
                    job_name_l.append(jname)
                    shell("""
                        if [ ! -d "{params.per_donor_log_dir}" ]; then mkdir -p {params.per_donor_log_dir}; fi
                        cmd_str="{i} {input.barcodes_vs_donor} {params.temp_bam_per_cell_dir} {params.split_bams_dir} {input.filt_bam} "
                        if [[ {params.filter_bed} != "None" && {params.filter_bed} != "False" ]]; then
                            if [[ {params.mito_info} != "None" && {params.mito_info} != "False" ]]; then
                                cmd_str+="{params.filter_bed} {params.mito_info} {params.chr_mito} "
                            else
                                cmd_str+="{params.filter_bed} "
                            fi
                        else
                            if [[ {params.mito_info} != "None" && {params.mito_info} != "False" ]]; then
                                cmd_str+="None {params.mito_info} {params.chr_mito} "
                            fi
                        fi
                        jid=$(bsub -J {j} -P acc_CommonMind -q express -n 1 -R span[hosts=1] \
                            -R rusage[mem=200] -W {params.time_limit_per_donor} \
                            -oo {params.per_donor_log_dir}{j}.stdout -eo {params.per_donor_log_dir}{j}.stderr \
                            -L /bin/bash "bash helper_sh_scripts/create_per_donor_bams.bash ${{cmd_str}} ")
                        jid=$(echo $jid | head -n1 | cut -d '<' -f2 | cut -d '>' -f1)
                        echo "Submitted script for donor {i} with jobid: ${{jid}}" >> {output[0]}
                        sleep 10
                    """, i=donor, j=jname)

            with open(output[0], "a") as fout:
                fout.write(f"Number of 'new' bam files expected at the completion of all the "
                    f"scripts for the bam file {input.filt_bam} is {len(job_name_l)}\n"
                )

        else:
            with open(output[0], "a") as fout:
                fout.write(f"Skipped bam file {input.filt_bam} as all "
                f"{len(hash_file.iloc[:, 0].unique())} donor file(s) was(were) "
                f"already present in the given output_folder {params.split_bams_dir}\n"
                )

