
rule Picard_GC_bias_metrics:
    input:
        bams=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}",
        bais=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bai']}"

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['gc_bias_metrics']}",
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['gc_summary_metrics']}"

    params:
        output_pref=lambda wildcards, output: output[0].replace(".txt", '.pdf'),
        window_size=config['picard_pipeline']["window_size"],
        genome_fasta=config['STARsolo_pipeline']['genome_pick']['fasta'] # config["genome_fasta"]

    resources:
        **get_resources("Picard_GC_bias_metrics")
       
    conda: "../envs/picard.yaml"

    envmodules:
        "picard/2.22.3"
     
    shell:
        """
        java -jar $PICARD CollectGcBiasMetrics I={input.bams} \
            O={output[0]} \
            CHART={params.output_pref} \
            SCAN_WINDOW_SIZE={params.window_size} S={output[1]} \
            R={params.genome_fasta}
        """
      


rule Picard_RNAseq_metrics:
    input:
        bams=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}",
        bais=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bai']}"

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['rnaseq_metrics']}"

    params:
        flat_ref=config['STARsolo_pipeline']['genome_pick']['flat_ref'], # config['picard_pipeline']['flat_ref'],
        strand=config['picard_pipeline']['strand']

    resources:
        **get_resources("Picard_RNAseq_metrics")

    conda: "../envs/picard.yaml"

    envmodules:
        "picard/2.22.3"

    shell:
        """
        java -jar $PICARD CollectRnaSeqMetrics I={input.bams} \
            O={output} REF_FLAT={params.flat_ref} \
            STRAND={params.strand}
        """
