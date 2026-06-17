# import os, uuid
import pandas as pd
from lib.io import (
    get_inp_splitBam, 
    get_bam,
    get_bam_to_split
)
from lib.params import (
    get_mito, 
    subset_to_chr, 
    get_mito_file,
)
from functools import partial


# Resource Allocation ------------------
def allocate_mem_FCB(wildcards, attempt):
    return 1000+200*(attempt-1)


def allocate_time_FCB(wildcards, attempt):
    return 20+20*(attempt-1)


def allocate_mem_SB(wildcards, attempt):
    return 250+100*(attempt-1)

# --------------------------------------

# FIX THE INPUT: INPUT SHOULD TAKE IN EITHER VIREO OR CALICO_SOLO'S H5AD AS INPUT, AS REQUIRED.
checkpoint create_inp_splitBams:
    input:
        partial(get_inp_splitBam, config=config)

    # priority: 7

    output:
        (
            f"{config['split_bams_pipeline']['inp_split_bams_dir']}"
            f"{config['fold_struct_bam_split1']}_bc_hash.txt"
        )

    params:
        conv=config['split_bams_pipeline']['donor_name_converter']['file'],
        from_col=config['split_bams_pipeline']['donor_name_converter']['from_column'],
        to_col=config['split_bams_pipeline']['donor_name_converter']['to_column'],
        demux_method=config['split_bams_pipeline']['split_by']['demux'],
        inp_ext=config['split_bams_pipeline']['split_by']['input'],
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
        unpack(partial(get_bam, config=config)),
        hash_file=(
            f"{config['split_bams_pipeline']['inp_split_bams_dir']}"
            f"{config['fold_struct_bam_split1']}_bc_hash.txt"
            )

    output:
        bam=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['filt_bam']}"
        ),
        bai=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['filt_bam']}.bai"
        )

    params:
        temp_bc=lambda wc: f"{config['split_bams_pipeline']['sort_temp_dir']}{{pool}}_bc.txt",
        threads=lambda _, resources: max(resources.cpus_per_task*6, 10)


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
        r"""
        mkdir -p {config[split_bams_pipeline][sort_temp_dir]}
        cut -f2 <(tail -n +2 {input.hash_file}) > {params.temp_bc}
        samtools view -@ {params.threads} -q 255 -D CB:{params.temp_bc} \
            {input.bam} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        rm {params.temp_bc}
        sleep 10
        """


rule filt_chr_bams:
    input:
        unpack(partial(get_bam, config=config))

    output:
        bam=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['short_bam']}"
        ), # generalize this
        bai=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['short_bam']}.bai"
        )

    params:
        sub_chr=partial(subset_to_chr, config=config),
        threads=lambda _, resources: max(resources.cpus_per_task*6, 10)

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
        samtools view -@ {params.threads} {input.bam} {params.sub_chr} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        """


rule filt_chr_bams_multiome:
    input:
        unpack(partial(get_bam, config=config))

    output:
        gex_bam=(
            f"{config['cellranger_arc_count']['bams_dir']}"
            f"{config['fold_struct']}{{bam}}"
            f"{config['split_bams_pipeline']['short_bam']}"
        ), # generalize this
        gex_bai=(
            f"{config['cellranger_arc_count']['bams_dir']}"
            f"{config['fold_struct']}{{bam}}"
            f"{config['split_bams_pipeline']['short_bam']}.bai"
        ),
        atac_bam=(
            f"{config['cellranger_arc_count']['bams_dir']}"
            f"{config['fold_struct']}{{bam}}"
            f"{config['split_bams_pipeline']['short_bam']}"
        ), # generalize this
        atac_bai=(
            f"{config['cellranger_arc_count']['bams_dir']}"
            f"{config['fold_struct']}{{bam}}"
            f"{config['split_bams_pipeline']['short_bam']}.bai"
        ) 

    params:
        sub_chr=partial(subset_to_chr, config=config),
        threads=lambda _, resources: max(resources.cpus_per_task*6, 10)

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
        samtools view -@ {params.threads} {input.bam} {params.sub_chr} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        """


rule create_bed:
    input:
        config['STARsolo_pipeline']['genome_pick']['fasta'] # config['genome_fasta']

    params:
        chr_prefix= lambda _: (
            '' if config['chr_prefix'] is None 
            else config['chr_prefix']
        )

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


def get_real_donors(wildcards):
    import pandas as pd
    from itertools import repeat

    ckpt = checkpoints.create_inp_splitBams.get(**wildcards)
    demux_out = ckpt.output[0]
    min_cells_threshold = config['split_bams_pipeline']['min_cells_per_donor']

    df = pd.read_csv(demux_out, sep="\t")

    # filter donors with meaningful data
    df = df.groupby('Subj_ID')['barcodes'].count()
    donors = df[df>min_cells_threshold].index.tolist()

    return {'pool': list(repeat(wildcards['pool'], len(donors))) , 'donor': donors}


# rule split_bams:
#     input:
#         filt_bam=get_bam_to_split,
#         barcodes_vs_donor=f"{config['split_bams_pipeline']['inp_split_bams_dir']}{config['fold_struct_bam_split1']}_bc_hash.txt" # with headers

#     params:
#         # split_at=config['split_bams_pipeline']['bc_per_donor'], # Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
#         # temp_dir=f"{config['split_bams_pipeline']['temp_dir']}",
#         # temp_bam=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}{{pool}}_{{donor}}.bam", # wildcards
#         temp_bam_per_cell_dir=f"{config['split_bams_pipeline']['new_temp_dir']}{config['fold_struct_bam_split2']}",
#         split_bams_dir=f"{config['split_bams_pipeline']['split_bams_dir']}{config['fold_struct_bam_split2']}",
#         per_donor_log_dir=config['split_bams_pipeline']['per_donor_split_log_dir'],
#         time_limit_per_donor=config['split_bams_pipeline']['time_per_donor'],
#         chr_mito=get_mito,
#         gt_check=config['gt_check'],
#         mito_info=get_mito_file,
#         filter_bed=None if config['gt_check'] else config['reg_chr_bed']

#     # For snakemake < v8
#     # threads: 1

#     resources:
#         cpus_per_task=1, # For snakemake > v8
#         mem_mb=allocate_mem_SB,
#         time_min=30

#     output:
#         # f"{config['split_bams_pipeline']['split_bams_proxy_dir']}{config['fold_struct_bam_split1']}.txt"  # Proxy to the output

#     run:
#         hash_file=pd.read_csv(input.barcodes_vs_donor, sep='\t')
#         proc_donors=[d for d in hash_file.iloc[:, 0].unique() if not os.path.isfile(os.path.join(params.split_bams_dir, f"{d}.bam.bai"))]
#         job_name_l=[]
#         samp_num=('-'.join(os.path.basename(input.filt_bam).replace('_filt.bam', '').split('-')[1:3])) # project-specific
#         # If proc_donors is not empty then run
#         if proc_donors:
#             with open(output[0], "a") as fout:
#                 fout.write("Going to process donor(s): {}\n".format(','.join(proc_donors)))

#             if params.gt_check:
#                 for donor in proc_donors:
#                     # jname=f"{samp_num}_" + list(shell("uuidgen", iterable=True))[0]
#                     jname=f"{samp_num}_{uuid.uuid4().hex}"
#                     job_name_l.append(jname)
#                     shell("""
#                         if [ ! -d "{params.per_donor_log_dir}" ]; then mkdir -p {params.per_donor_log_dir}; fi
#                         jid=$(bsub -J {j} -P acc_CommonMind -q express -n 1 -R span[hosts=1] \
#                             -R rusage[mem=200] -W {params.time_limit_per_donor} \
#                             -oo {params.per_donor_log_dir}{j}.stdout -eo {params.per_donor_log_dir}{j}.stderr \
#                             -L /bin/bash "bash helper_sh_scripts/create_per_donor_bams.bash {i} \
#                             {input.barcodes_vs_donor} {params.temp_bam_per_cell_dir} \
#                             {params.split_bams_dir} {input.filt_bam} ")
#                         jid=$(echo $jid | head -n1 | cut -d '<' -f2 | cut -d '>' -f1)
#                         echo "Submitted script for donor {i} with jobid: ${{jid}}" >> {output[0]}
#                         sleep 10
#                     """, i=donor, j=jname)
#             else:
#                 for donor in proc_donors:
#                     jname=f"{samp_num}_" + list(shell("uuidgen", iterable=True))[0]
#                     job_name_l.append(jname)
#                     shell("""
#                         if [ ! -d "{params.per_donor_log_dir}" ]; then mkdir -p {params.per_donor_log_dir}; fi
#                         cmd_str="{i} {input.barcodes_vs_donor} {params.temp_bam_per_cell_dir} {params.split_bams_dir} {input.filt_bam} "
#                         if [[ {params.filter_bed} != "None" && {params.filter_bed} != "False" ]]; then
#                             if [[ {params.mito_info} != "None" && {params.mito_info} != "False" ]]; then
#                                 cmd_str+="{params.filter_bed} {params.mito_info} {params.chr_mito} "
#                             else
#                                 cmd_str+="{params.filter_bed} "
#                             fi
#                         else
#                             if [[ {params.mito_info} != "None" && {params.mito_info} != "False" ]]; then
#                                 cmd_str+="None {params.mito_info} {params.chr_mito} "
#                             fi
#                         fi
#                         jid=$(bsub -J {j} -P acc_CommonMind -q express -n 1 -R span[hosts=1] \
#                             -R rusage[mem=200] -W {params.time_limit_per_donor} \
#                             -oo {params.per_donor_log_dir}{j}.stdout -eo {params.per_donor_log_dir}{j}.stderr \
#                             -L /bin/bash "bash helper_sh_scripts/create_per_donor_bams.bash ${{cmd_str}} ")
#                         jid=$(echo $jid | head -n1 | cut -d '<' -f2 | cut -d '>' -f1)
#                         echo "Submitted script for donor {i} with jobid: ${{jid}}" >> {output[0]}
#                         sleep 10
#                     """, i=donor, j=jname)

#             with open(output[0], "a") as fout:
#                 fout.write(f"Number of 'new' bam files expected at the completion of all the "
#                     f"scripts for the bam file {input.filt_bam} is {len(job_name_l)}\n"
#                 )

#         else:
#             with open(output[0], "a") as fout:
#                 fout.write(f"Skipped bam file {input.filt_bam} as all "
#                 f"{len(hash_file.iloc[:, 0].unique())} donor file(s) was(were) "
#                 f"already present in the given output_folder {params.split_bams_dir}\n"
#                 )


# ---------------------------------------------------------------------------
# Rule: one job per donor, both gt_check modes
# ---------------------------------------------------------------------------

rule split_bams:
    input:
        unpack(partial(get_bam_to_split, config=config)),
        barcodes_vs_donor = (
            f"{config['split_bams_pipeline']['inp_split_bams_dir']}"
            f"{config['fold_struct_bam_split1']}_bc_hash.txt"
        )
    output:
        bam = (
            f"{config['split_bams_pipeline']['split_bams_dir']}"
            f"{config['fold_struct_bam_split2']}"
            "{donor}.bam"
        ),
        bai = (
            f"{config['split_bams_pipeline']['split_bams_dir']}"
            f"{config['fold_struct_bam_split2']}"
            "{donor}.bam.bai"
        )
    params:
        temp_dir      = (
            f"{config['split_bams_pipeline']['new_temp_dir']}"
            f"{config['fold_struct_bam_split2']}"
        ),
        gt_check      = config["gt_check"],
        chr_mito      = partial(get_mito, config=config),
        mito_info     = partial(get_mito_file, config=config),
        filter_bed    = lambda _: (
            None if config["gt_check"]
            else config.get("reg_chr_bed")
        ),
        script        = "helper_sh_scripts/create_per_donor_bams.bash"
    resources:
        cpus_per_task = 1,
        mem_mb        = allocate_mem_SB,
        time_min      = 30
    log:
        (
            f"{config["split_bams_pipeline"]["per_donor_split_log_dir"]}"
            f"{config['fold_struct_bam_split2']}"
            "{donor}.log",
        )
    shell:
        """
        # Build the optional-argument string based on mode.
        # gt_check mode: no BED or mito refinement needed.
        # full-split mode: pass BED and/or mito file if configured.

        EXTRA_ARGS=""

        if [[ "{params.gt_check}" != "True" ]]; then
            if [[ "{params.filter_bed}" != "None" && "{params.filter_bed}" != "False" ]]; then
                EXTRA_ARGS+=" --bed {params.filter_bed}"
            fi
            if [[ "{params.mito_info}" != "None" && "{params.mito_info}" != "False" ]]; then
                EXTRA_ARGS+=" --mito_file {params.mito_info} --mito_pref {params.chr_mito}"
            fi
        fi

        bash {params.script} \
            --donor        {wildcards.donor}         \
            --hash_file    {input.barcodes_vs_donor} \
            --temp_dir     {params.temp_dir}         \
            --out_dir      $(dirname {output.bam})   \
            --bam          {input.bam}          \
            ${{EXTRA_ARGS}}                          \
            > {log} 2>&1
        """
