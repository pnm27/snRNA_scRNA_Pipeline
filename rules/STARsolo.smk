import os, glob2, re


def get_r1_fastqs(wildcards):
    temp_fold=f"{config['fold_struct']}".format(**wildcards) # Add wildcards
    r1_files = sorted(glob2.glob(f"{config['cDNA_fastqs_dir']}{temp_fold}*{config['R1_suffix']}"))
    return r1_files


def get_r2_fastqs(wildcards):
    temp_fold=f"{config['fold_struct']}".format(**wildcards)
    r2_files = sorted(glob2.glob(f"{config['cDNA_fastqs_dir']}{temp_fold}*{config['R2_suffix']}"))
    return r2_files


# To check if the last digit of the line in the error log of STARsolo is a number
def check_isnumber(x):
    try:
        int(x)
        return True
   
    except ValueError:
        return False


# We can use this to similarly change other params in the log
def get_limitsjdbval_coll(wildcards, resources):
    '''
    This function reads the log file created per attempt to change the parameter "limitsjdbInsertNsj" and "limitOutSJcollapsed" in STARsolo
    Serially produce output as a list in the sequence "limitsjdbInsertNsj", "limitOutSJcollapsed", etc.
    '''
    # This is to check the log file produced after each attempt for the error value
    file_p_temp = f"{config['fold_struct']}".format(**wildcards)
    log_list = glob2.glob(
        "{}{}_STARsolo_log.txt*".format(config['STARsolo_pipeline']['bams_dir'], 
            file_p_temp)
        )
    ins_nsj = 1000000 # DEFAULT
    sj_collap = 1000000 # DEFAULT
    limitbamsortram = 5000000000
    # limitbamsortram = int(resources.mem_mb * 0.5 * 1000000) # resources.mem_mb * (resources.cpus_per_task - 1) * 1000000 # DEFAULT
    for log_file in log_list: 
        with open(log_file) as fin:
            for line in fin:
                if line.lower().startswith("solution") and "limitSjdbInsertNsj" in line \
                    and check_isnumber(line.split()[-1]):
                    # print("Found an Error with the parameter limitSjdbInsertNsj. Changing defaulkt values of the parameters \
                        # \"limitSjdbInsertNsj\" and \"limitOutSJcollapsed\" from the default value of 1000000 to {}".format(line.split()[-1]))
                    ins_nsj = int(line.split()[-1]) if int(line.split()[-1]) > ins_nsj else ins_nsj
                    sj_collap = ins_nsj

                elif line.lower().startswith("solution") and "limitOutSJcollapsed" in line:
                    # print("Found an Error with limitOutSJcollapsed. Changing from the default value of 1000000 to {}".format(1000000*(1+resources.attempt)))
                    sj_collap = 1000000*(1+resources.attempt)
                    ins_nsj = sj_collap
                
                elif line.lower().startswith("solution") and "limitBAMsortRAM" in line:
                    # print("Found an Error with limitOutSJcollapsed. Changing from the default value of 1000000 to {}".format(1000000*(1+resources.attempt)))
                    val = re.search(r"--limitBAMsortRAM ([0-9]+) ", line)
                    limitbamsortram = int(val.group(1)) if val is not None else limitbamsortram

                else:
                    continue

            # Empty but existing file
            else:
                continue
                    

   # This is to check the parameters file, else block WILL execute after the for block as the for block as no break in it (need revision)
   # Use try except block to catch file issues
    else:
        if os.path.isfile("{}{pool}-cDNA.txt".format(config['STARsolo_pipeline']['star_params_dir'], **wildcards)): # wildcard
            with open("{}{pool}-cDNA.txt".format(config['STARsolo_pipeline']['star_params_dir'], **wildcards)) as fin: # wildcard
                for line in fin:
                    # print("Found values of \"limitSjdbInsertNsj\" and \"limitOutSJcollapsed\" from the previous successfull run in {}. Using the same value".format(config['star_params_dir']))
                    val = re.search(r"--limitSjdbInsertNsj ([0-9]+) ", line)
                    temp_nsj = int(val.group(1)) if val is not None else 0
                    val = re.search(r"--limitOutSJcollapsed ([0-9]+) ", line)
                    temp_sj_coll = int(val.group(1)) if val is not None else 0
                    ins_nsj = max(temp_nsj, temp_sj_coll, ins_nsj, sj_collap)
                    sj_collap = ins_nsj
                    val = re.search(r"--limitBAMsortRAM ([0-9]+) ", line)
                    limitbamsortram = int(val.group(1)) if val is not None else 0


    return [ins_nsj, sj_collap, limitbamsortram]


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
        R1=get_r1_fastqs,
        R2=get_r2_fastqs

    # priority: 10

    params:
        # gtf=config['gtf_file'],
        gtf=config['STARsolo_pipeline']['genome_pick']['gtf'], # USED to be ['gtf_file']
        genome_dir=config['STARsolo_pipeline']['genome_pick']['genome'], # USED to be ['STARsolo_pipeline']['genome_dir']
        overhang=config['STARsolo_pipeline']['genome_pick']['overhang'], # USED to be ['STARsolo_pipeline']['sjdboverhang']
        opt_params=get_limitsjdbval_coll,
        chemistry=config['STARsolo_pipeline']['soloType'], # For STARsolo
        whitelist=config['whitelist'], # V3 whitelist
        UMI_length=config['STARsolo_pipeline']['umi_len'], # V3 
        SAM_attr=config['STARsolo_pipeline']['SAM_attr'],
        features=config['STARsolo_pipeline']['features'],
        save_params=f"{config['STARsolo_pipeline']['star_params_dir']}{{pool}}-cDNA.txt",  # WILDCARDS
        star_def_log_out=lambda wildcards, output: output.bam.replace(config['STARsolo_pipeline']['bam'], "_Log.out"),
        solo_cell_filter=config['STARsolo_pipeline']['solo_cell_filter'],
        out_pref=lambda wildcards, output: output.bam.replace(config['STARsolo_pipeline']['bam'], '_'),
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
        r1=$(echo "{input.R1}" | tr '[:blank:]' ',')
        r2=$(echo "{input.R2}" | tr '[:blank:]' ',')
        echo "{params.opt_params[0]}, {params.opt_params[1]}, {params.opt_params[2]}, {resources.attempt}"
        if [ ! -d {config[STARsolo_pipeline][star_params_dir]} ]; then mkdir -p {config[STARsolo_pipeline][star_params_dir]}; fi
        if [[ "{config[STARsolo_pipeline][extra_params]}" == "None" ]]; then
            STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[0]} \
            --twopassMode Basic --readFilesCommand zcat --readFilesIn ${{r2}} ${{r1}} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} \
            --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} \
            --limitOutSJcollapsed {params.opt_params[1]} --outSAMtype BAM SortedByCoordinate --runThreadN {params.threads} --outFileNamePrefix {params.out_pref} \
            --limitBAMsortRAM {params.opt_params[2]} --outBAMsortingBinsN 50 &> {log}_{resources.attempt}
        else
            STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[0]} \
            --twopassMode Basic --readFilesCommand zcat --readFilesIn ${{r2}} ${{r1}} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} \
            --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} \
            --limitOutSJcollapsed {params.opt_params[1]} --outSAMtype BAM SortedByCoordinate --runThreadN {params.threads} --outFileNamePrefix {params.out_pref} \
            --limitBAMsortRAM {params.opt_params[2]} --outBAMsortingBinsN 50 {config[STARsolo_pipeline][extra_params]} &> {log}_{resources.attempt}
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
