rule Picard_GC_bias_metrics:
    input:
        bams=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"

    # priority: 9

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['gc_bias_metrics']}",
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['gc_summary_metrics']}"

    params:
        output_pref=lambda wildcards, output: output[0].replace(f"{config['picard_pipeline']['gc_bias_metrics']}", '_'),
        window_size=config['picard_pipeline']["window_size"],
        genome_fasta=config["genome_fasta"]

    resources:
        mem_mb=allocate_mem_PGCB,
        time_min=allocate_time_PGCB
       
    threads: 2

    conda: "../envs/picard.yaml"

    envmodules:
        "picard/2.22.3"
     
    shell:
        """
        ml picard/2.22.3
        java -jar $PICARD CollectGcBiasMetrics I={input.bams} O={output[0]} CHART={params.output_pref}gc_bias_metrics.pdf SCAN_WINDOW_SIZE={params.window_size} S={output[1]} R={params.genome_fasta}
        """
      


rule Picard_RNAseq_metrics:
    input:
        bams=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['bam']}"

    # priority: 9

    output:
        f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['picard_pipeline']['rnaseq_metrics']}"

    params:
        flat_ref=config['picard_pipeline']['flat_ref'],
        strand=config['picard_pipeline']['strand']

    resources:
        mem_mb=allocate_mem_PRNA,
        time_min=allocate_time_PRNA

    threads: 4

    conda: "../envs/picard.yaml"

    envmodules:
        "picard/2.22.3"

    shell:
        """
        ml picard/2.22.3
        java -jar $PICARD CollectRnaSeqMetrics I={input.bams} O={output} REF_FLAT={params.flat_ref} STRAND={params.strand}
        """
