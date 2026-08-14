
from lib.io import get_r1_fastqs, get_r2_fastqs
from functools import partial
from lib.params import get_limitsjdbval_coll



# Resource Allocation ------------------

def allocate_mem_SS(wildcards, attempt):
    if not SS_MEM or wildcards.pool not in SS_MEM:
        return 45000+1000*(attempt-1) # DEFAULT
    else:
        return SS_MEM[wildcards.pool] +4000*(attempt-1)

def allocate_time_SS(wildcards, attempt):
    return 1440

# --------------------------------------


rule STARsolo_sort:
    input:
        R1=partial(get_r1_fastqs, config=config),
        R2=partial(get_r2_fastqs, config=config)

    # priority: 10

    params:
        # gtf=config['gtf_file'],
        ss_style_r1=lambda wildcards, input: ",".join(input.R1),
        ss_style_r2=lambda wildcards, input: ",".join(input.R2),
        gzip_inputs=lambda wildcards, input: ( 
            "--readFilesCommand zcat" 
            if input.R1[0].endswith(".gz") 
            else ""
        ),
        gtf=config['STARsolo_pipeline']['genome_pick']['gtf'], # USED to be ['gtf_file']
        genome_dir=config['STARsolo_pipeline']['genome_pick']['genome'], # USED to be ['STARsolo_pipeline']['genome_dir']
        overhang=config['STARsolo_pipeline']['genome_pick']['overhang'], # USED to be ['STARsolo_pipeline']['sjdboverhang']
        opt_params=lambda wildcards, resources: (
            partial(get_limitsjdbval_coll, 
                config=config, resources=resources)
        ), # To change STARsolo params based on the log file of the previous attempt
        chemistry=config['STARsolo_pipeline']['soloType'], # For STARsolo
        whitelist=config['whitelist'], # V3 whitelist
        UMI_length=config['STARsolo_pipeline']['umi_len'], # V3 
        SAM_attr=config['STARsolo_pipeline']['SAM_attr'],
        features=config['STARsolo_pipeline']['features'],
        save_params=f"{config['STARsolo_pipeline']['star_params_dir']}{{pool}}-cDNA.txt",  # WILDCARDS
        star_def_log_out=lambda wildcards, output: (
            output.bam.replace(
                config['STARsolo_pipeline']['bam'], 
                "_Log.out"
            )
        ),
        solo_cell_filter=config['STARsolo_pipeline']['solo_cell_filter'],
        out_pref=lambda wildcards, output: (
            output.bam.replace(
                config['STARsolo_pipeline']['bam'], 
                '_'
            )
        ),
        threads=config['STARsolo_pipeline']['run_threads']


    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['STAR_log_final']}",
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_summary']}",
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['gene_features']}",
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['gene_summary']}",
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['barcodes_stats']}",
        gf_feat=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_features']}",
        bam=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}",
        gf_mat=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_matrix']}",
        gf_bc=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_barcodes']}"

    log:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}_STARsolo_log.txt"
   
    resources:
        cpus_per_task=10, # For snakemake > v8
        mem_mb=allocate_mem_SS,
        time_min=allocate_time_SS,
        attempt=lambda wildcards, attempt: attempt

    # For snakemake < v8
    # threads: 6

    envmodules:
        f"{config['STAR_version']}"

    shell:
        r"""
        echo "Attempt: {resources.attempt}"
        if [ ! -d {config[STARsolo_pipeline][star_params_dir]} ]; then mkdir -p {config[STARsolo_pipeline][star_params_dir]}; fi
        if [[ "{config[STARsolo_pipeline][extra_params]}" == "None" ]]; then
            STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[ins_nsj]} \
            --twopassMode Basic {params.gzip_inputs} --readFilesIn {params.ss_style_r2} {params.ss_style_r1} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} \
            --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} \
            --limitOutSJcollapsed {params.opt_params[sj_collap]} --outSAMtype BAM SortedByCoordinate --runThreadN {params.threads} --outFileNamePrefix {params.out_pref} \
            --limitBAMsortRAM {params.opt_params[limitbamsortram]} --outBAMsortingBinsN 50 &> {log}_{resources.attempt}
        else
            STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[ins_nsj]} \
            --twopassMode Basic {params.gzip_inputs} --readFilesIn {params.ss_style_r2} {params.ss_style_r1} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} \
            --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} \
            --limitOutSJcollapsed {params.opt_params[sj_collap]} --outSAMtype BAM SortedByCoordinate --runThreadN {params.threads} --outFileNamePrefix {params.out_pref} \
            --limitBAMsortRAM {params.opt_params[limitbamsortram]} --outBAMsortingBinsN 50 {config[STARsolo_pipeline][extra_params]} &> {log}_{resources.attempt}
        fi
        files=( {output.gf_feat} {output.gf_mat} {output.gf_bc} )
        for i in ${{files[@]}}; do
            base=${{i%".gz"}}
            if [ -f "${{base}}" ] && [ ! -f "${{i}}" ]; then
                gzip "${{base}}"
            fi
        done
        a=$(grep -n "^##### Final effective command line" {params.star_def_log_out} | cut -d ":" -f1)
        a=$((a+1))
        if [ ! -f "{params.save_params}" ]; then 
            sed -n "${{a}}p" {params.star_def_log_out} > {params.save_params}
        else 
            sed -n "${{a}}p" {params.star_def_log_out} > {params.save_params}_{resources.attempt}
            if cmp --silent {params.save_params} {params.save_params}_{resources.attempt}; then
                rm {params.save_params}_{resources.attempt}
            else
                rm {params.save_params}
                mv {params.save_params}_{resources.attempt} {params.save_params}
            fi
        fi
        op_dir={output.bam}
        op_dir=$(dirname ${{op_dir}})
        find ${{op_dir}} -type f -name '*mate*' -print0 | \
        tar --null --remove-files -zcvf {params.out_pref}Unmapped_fastqs.tar.gz --files-from=-
        """


rule samtools_index:
    input:
        bam=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"

    params:
        threads=lambda wildcards, resources: resources.cpus_per_task*4-2

    output:
        bai=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bai']}"

    resources:
        cpus_per_task=4, # For snakemake > v8
        mem_mb=5000,
        time_min=720,
        attempt=lambda wildcards, attempt: attempt
    
    log:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}_samtools_index_log.txt"

    envmodules:
        "samtools/1.21"

    shell:
        """
        samtools index -@ {params.threads} {input.bam} &>> {log}
        """
