# Vireo output has a fixed 37 comments with one header line i.e. 38 lines from the top are of no analytical use
rule calc_vireoSNP_vals:
    input:
        f"{config['phe_demux_pipeline']['vireosnp_dir']}{config['fold_struct_phe_demux']}{{vcf_type}}/{config['phe_demux_pipeline']['donors_vcf']}",
        f"{config['phe_demux_pipeline']['vireosnp_dir']}{config['fold_struct_phe_demux']}{{vcf_type}}/{config['phe_demux_pipeline']['donors_classification']}"
  
    output:
        f"{config['phe_demux_pipeline']['vireosnp_dir']}{config['fold_struct_phe_demux']}{{vcf_type}}/{config['phe_demux_pipeline']['donors_vcf']}"
        #f"{config['filt_vcf_dir']}{config['fold_struct_kb']}{config['filt_vcf']}"


    # group: "phenotype-demux"

    # params:
    #     donor_info=get_donors,
    #     geno_tag=config['phe_demux_pipeline']['donor_genotype'],
    #     output_prefix=lambda wildcards, output: output[0].replace(f"/{config['phe_demux_pipeline']['donors_vcf']}", '')

    threads: 6

    resources:
        mem_mb=allocate_mem_cvv,
        time_min=allocate_time_cvv
        
    shell:
        """
        if [ "{params.donor_info[donor]}" -eq "6" ]; then
            vireo -c {input} -d {params.donor_info[vcf]} -N 6 -o {params.output_prefix} -t {params.geno_tag} --noPlot -p 20 --forceLearnGT
        else
            vireo -c {input} -d {params.donor_info[vcf]} -N 6 -o {params.output_prefix} -t {params.geno_tag} --noPlot -p 20
        fi
        sleep 180
        """