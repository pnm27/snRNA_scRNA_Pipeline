
from lib.io import get_r1_fastqs, get_r2_fastqs, pool_vcf, wasp_vcf
from functools import partial
from lib.params import get_limitsjdbval_coll
from lib.aligner import starsolo_matrix_keys



WASP = config['modules']['STARsolo']['wasp']


rule STARsolo_sort:
    input:
        R1=partial(get_r1_fastqs, config=config),
        R2=partial(get_r2_fastqs, config=config),
        # WASP mode: the pool's vcf (rule 'wasp_vcf')
        wasp_vcf=wasp_vcf(config) if WASP else []

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
        gtf=config['modules']['STARsolo']['genome_pick']['gtf'], # USED to be ['gtf_file']
        genome_dir=config['modules']['STARsolo']['genome_pick']['genome'], # USED to be ['STARsolo_pipeline']['genome_dir']
        overhang=config['modules']['STARsolo']['genome_pick']['overhang'], # USED to be ['STARsolo_pipeline']['sjdboverhang']
        opt_params=lambda wildcards, resources: (
            partial(get_limitsjdbval_coll, 
                config=config, resources=resources)
        ), # To change STARsolo params based on the log file of the previous attempt
        chemistry=config['modules']['STARsolo']['soloType'], # For STARsolo
        whitelist=config['references']['whitelist'], # V3 whitelist
        UMI_length=config['modules']['STARsolo']['umi_len'], # V3 
        SAM_attr=(
            config['modules']['STARsolo']['SAM_attr']
            + (" vW" if WASP and "vW" not in config['modules']['STARsolo']['SAM_attr'].split() else "")
        ),
        wasp=lambda wildcards, input: (
            f"--waspOutputMode SAMtag --varVCFfile {input.wasp_vcf}" if WASP else ""
        ),
        features=config['modules']['STARsolo']['features'],
        save_params=f"{config['modules']['STARsolo']['star_params_dir']}{{pool}}-cDNA.txt",  # WILDCARDS
        star_def_log_out=lambda wildcards, output: (
            output.bam.replace(
                config['modules']['STARsolo']['bam'], 
                "_Log.out"
            )
        ),
        solo_cell_filter=config['modules']['STARsolo']['solo_cell_filter'],
        out_pref=lambda wildcards, output: (
            output.bam.replace(
                config['modules']['STARsolo']['bam'], 
                '_'
            )
        ),
        threads=config['modules']['STARsolo']['run_threads'],
        gene_files=lambda wildcards, output: " ".join(
            output[k] for k in ['gene_mat', 'gene_bc'] if k in output.keys()
        )


    output:
        f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['STAR_log_final']}",
        f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['genefull_summary']}",
        f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['gene_summary']}",
        f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['barcodes_stats']}",
        g_feat=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['gene_features']}",
        gf_feat=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['genefull_features']}",
        bam=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['bam']}",
        gf_mat=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['genefull_matrix']}",
        gf_bc=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['genefull_barcodes']}",
        # Gene's matrix and barcodes, when it's the count matrix
        **{
            name: f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo'][key]}"
            for name, key in zip(['gene_mat', 'gene_bc'], starsolo_matrix_keys(config)[:2])
            if key.startswith('gene_')
        }

    log:
        f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}_STARsolo_log.txt"
   
    resources:
        **get_resources("STARsolo_sort"),
        attempt=lambda wildcards, attempt: attempt

    envmodules:
        f"{config['references']['tools']['star']}"

    shell:
        r"""
        echo "Attempt: {resources.attempt}"
        if [ ! -d {config[modules][STARsolo][star_params_dir]} ]; then mkdir -p {config[modules][STARsolo][star_params_dir]}; fi
        if [[ "{config[modules][STARsolo][extra_params]}" == "None" ]]; then
            STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[ins_nsj]} \
            --twopassMode Basic {params.gzip_inputs} --readFilesIn {params.ss_style_r2} {params.ss_style_r1} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} \
            --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} \
            --limitOutSJcollapsed {params.opt_params[sj_collap]} --outSAMtype BAM SortedByCoordinate --runThreadN {params.threads} --outFileNamePrefix {params.out_pref} \
            --limitBAMsortRAM {params.opt_params[limitbamsortram]} --outBAMsortingBinsN 50 {params.wasp} &> {log}_{resources.attempt}
        else
            STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[ins_nsj]} \
            --twopassMode Basic {params.gzip_inputs} --readFilesIn {params.ss_style_r2} {params.ss_style_r1} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} \
            --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} \
            --limitOutSJcollapsed {params.opt_params[sj_collap]} --outSAMtype BAM SortedByCoordinate --runThreadN {params.threads} --outFileNamePrefix {params.out_pref} \
            --limitBAMsortRAM {params.opt_params[limitbamsortram]} --outBAMsortingBinsN 50 {params.wasp} {config[modules][STARsolo][extra_params]} &> {log}_{resources.attempt}
        fi
        # gzip the matrix files (paths ending in .gz in the config)
        files=( {output.gf_feat} {output.gf_mat} {output.gf_bc} {output.g_feat} {params.gene_files} )
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
        bam=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['bam']}"

    params:
        threads=lambda wildcards, resources: resources.cpus_per_task*4-2

    output:
        bai=f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}{config['modules']['STARsolo']['bai']}"

    resources:
        **get_resources("samtools_index"),
        attempt=lambda wildcards, attempt: attempt
    
    log:
        f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}_samtools_index_log.txt"

    envmodules:
        "samtools/1.21"

    shell:
        """
        samtools index -@ {params.threads} {input.bam} &>> {log}
        """


# STAR reads the genotype of the vcf's first sample only, so for a pool of
# donors every biallelic SNV is written as heterozygous (0|1): WASP then
# checks the reads over all the pool's variants
rule wasp_vcf:
    input:
        partial(pool_vcf, config=config)

    output:
        temp(wasp_vcf(config))

    resources:
        **get_resources("wasp_vcf")

    shell:
        r"""
        zcat -f {input} | awk 'BEGIN {{OFS="\t"}}
            /^##/ {{ next }}
            /^#/ {{ print "##fileformat=VCFv4.2"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","POOL"; next }}
            length($4) == 1 && length($5) == 1 && $5 != "." {{ print $1,$2,$3,$4,$5,".",".",".","GT","0|1" }}' \
            > {output}
        """
