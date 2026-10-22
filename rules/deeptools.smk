import glob2, os


def get_BC_input(wildcards):
    # WILDCARDS
    all_bams = sorted(glob2.glob(
        f"{config['modules']['deeptools']['multiBamSummary']['bams_dir']}{config['modules']['deeptools']['fold_struct']}*{config['experiment']['inputs']['R1_suffix']}/".format(id2=wildcards.id2) + 
        f"{config['modules']['deeptools']['fold_struct']}*{config['modules']['deeptools']['multiBamSummary']['bam']}".format(id2=wildcards.id2)))

    return all_bams


rule multiBamSummary:
    input:
        get_BC_input
    # priority: 9

    output:
        f"{config['modules']['deeptools']['multiBamSummary']['outdir']}{config['modules']['deeptools']['fold_struct']}.npz"
        
    # params:
    #     outfile_fmt=config['modules']['deeptools']['multiBamSummary']['out_format'] if not config['modules']['deeptools']['multiBamSummary']['out_format'].startswith('.') else config['modules']['deeptools']['multiBamSummary']['out_format'][1:]

    conda: "../envs/deeptools.yaml"

    resources:
        **get_resources("multiBamSummary")
       
    # group: "PICARD_metrics"
     
    shell:
        """
        multiBamSummary bins --bamfiles {input} -o {output} -v -p max --smartLabels
        """


rule plotCorrelation:
    input:
        f"{config['modules']['deeptools']['multiBamSummary']['outdir']}{config['modules']['deeptools']['fold_struct']}.npz"
    # priority: 9

    output:
        f"{config['modules']['deeptools']['plotCorrelation']['outdir']}{config['modules']['deeptools']['fold_struct']}{config['modules']['deeptools']['plotCorrelation']['out_fmt']}"
        
    params:
        correlation_method=config['modules']['deeptools']['plotCorrelation']['corr_method'],
        plot_type=config['modules']['deeptools']['plotCorrelation']['plot_type'],
        corr_matrix=lambda wildcards, output: output[0].replace(config['modules']['deeptools']['plotCorrelation']['out_fmt'], '.txt'),
        labels=lambda wildcards, input: [os.path.basename(i).replace('.bam', '') for i in input ]

    conda: "../envs/deeptools.yaml"

    resources:
        **get_resources("plotCorrelation")
       
    # group: "PICARD_metrics"
     
    shell:
        """
        plotCorrelation -in {input} -c {params.correlation_method} -p {params.plot_type} -o {output} --outFileCorMatrix {params.corr_matrix}
        """
        