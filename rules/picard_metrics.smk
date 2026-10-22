from lib.aligner import bam_paths, out_prefix


rule Picard_GC_bias_metrics:
    input:
        bams=bam_paths(config, 'cDNA')['bam'],
        bais=bam_paths(config, 'cDNA')['bai']

    output:
        f"{out_prefix(config)}{config['modules']['picard']['gc_bias_metrics']}",
        f"{out_prefix(config)}{config['modules']['picard']['gc_summary_metrics']}"

    params:
        output_pref=lambda wildcards, output: output[0].replace(".txt", '.pdf'),
        window_size=config['modules']['picard']["window_size"],
        genome_fasta=config['modules']['STARsolo']['genome_pick']['fasta'] # config["genome_fasta"]

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
        bams=bam_paths(config, 'cDNA')['bam'],
        bais=bam_paths(config, 'cDNA')['bai']

    output:
        f"{out_prefix(config)}{config['modules']['picard']['rnaseq_metrics']}"

    params:
        flat_ref=config['modules']['STARsolo']['genome_pick']['flat_ref'], # config['modules']['picard']['flat_ref'],
        strand=config['modules']['picard']['strand']

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
