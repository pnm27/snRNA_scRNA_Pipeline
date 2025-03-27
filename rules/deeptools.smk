import glob2, os


def get_BC_input(wildcards):
    # WILDCARDS
    all_bams = sorted(glob2.glob(
        f"{config['multiBamSummary']['bams_dir']}{config['fold_struct_deeptools']}*{config['R1_suffix']}/".format(id2=wildcards.id2) + \ 
        f"{config['fold_struct_deeptools']}*{config['multiBamSummary']['bam']}".format(id2=wildcards.id2)))

    return all_bams


# Resource Allocation ------------------
def allocate_mem_MBS(wildcards, attempt):
    return 400 + 200*(attempt-1)


def allocate_time_MBS(wildcards, attempt):
    return 300 + 30*(attempt-1)


def allocate_mem_PC(wildcards, attempt):
    return 200 + 500*(attempt-1)


def allocate_time_PC(wildcards, attempt):
    return 10 + 200*(attempt-1)

# --------------------------------------

rule multiBamSummary:
    input:
        get_BC_input
    # priority: 9

    output:
        f"{config['multiBamSummary']['outdir']}{config['fold_struct_deeptools']}.npz"
        
    # params:
    #     outfile_fmt=config['multiBamSummary']['out_format'] if not config['multiBamSummary']['out_format'].startswith('.') else config['multiBamSummary']['out_format'][1:]

    conda: "../envs/deeptools.yaml"

    resources:
        cpus_per_task=5, # For snakemake > v8
        mem_mb=allocate_mem_MBS,
        time_min=allocate_time_MBS
       
    # For snakemake < v8
    # threads: 5

    # group: "PICARD_metrics"
     
    shell:
        """
        multiBamSummary bins --bamfiles {input} -o {output} -v -p max --smartLabels
        """


rule plotCorrelation:
    input:
        f"{config['multiBamSummary']['outdir']}{config['fold_struct_deeptools']}.npz"
    # priority: 9

    output:
        f"{config['plotCorrelation']['outdir']}{config['fold_struct_deeptools']}{config['plotCorrelation']['out_fmt']}"
        
    params:
        correlation_method=config['plotCorrelation']['corr_method'],
        plot_type=config['plotCorrelation']['plot_type'],
        corr_matrix=lambda wildcards, output: output[0].replace(config['plotCorrelation']['out_fmt'], '.txt'),
        labels=lambda wildcards, input: [os.path.basename(i).replace('.bam', '') for i in input ]

    conda: "../envs/deeptools.yaml"

    resources:
        cpus_per_task=2, # For snakemake > v8
        mem_mb=allocate_mem_PC,
        time_min=allocate_time_PC
       
    # For snakemake < v8
    # threads: 2

    # group: "PICARD_metrics"
     
    shell:
        """
        plotCorrelation -in {input} -c {params.correlation_method} -p {params.plot_type} -o {output} --outFileCorMatrix {params.corr_matrix}
        """