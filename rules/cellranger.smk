import os, glob2, re


def get_fastqs(wildcards):
    temp_fold=f"{config['fold_struct']}".format(**wildcards) # Add wildcards
    all_files = []
    all_files.extend(glob2.glob(f"{config['cDNA_fastqs_dir']}"
                f"{temp_fold}*{config['multiome_suffix']}"))
    all_files.extend(glob2.glob(f"{config['ATAC_fastqs_dir']}"
                f"{temp_fold}*{config['multiome_suffix']}"))

    return sorted(all_files)


# Resource Allocation ------------------

def allocate_mem_ICAC(wildcards, attempt):
    return 2000+1000*(attempt-1)

def allocate_time_ICAC(wildcards, attempt):
    return 5

def allocate_mem_CAC(wildcards, attempt):
    return 7000+2000*(attempt-1)

def allocate_time_CAC(wildcards, attempt):
    return 2880

# --------------------------------------

rule inp_cellranger_arc_count:
    input:
        get_fastqs

    params:
        samp_id=lambda wildcards: wildcards.pool

    output:
        "{pool}_metadata.csv"

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_ICAC,
        time_min=allocate_time_ICAC

    conda: "../envs/basic_sctools.yaml"

    shell:
        """
         python3 helper_py_scripts/create_inp_cellrangerArcCount.py \
            --output {output} -s {params.samp_id} {input} 
        """


rule cellranger_arc_count:
    """
    Whole directory is added as output to ensure all files are protected.
    This is just a sloppy way of not listing all cellranger-arc count outputs.
    This also ensures that, once write protection is removed, the whole directory 
    is wiped out before trying to re-run cellranger-arc count.

    Pipestance contains useful information.
    """
    input:
        "{pool}_metadata.csv",
        get_fastqs
    # priority: 10

    params:
        ref=config['cellranger_arc_ref'],
        max_localcores=lambda wildcards, resources: resources.cpus_per_task*4,
        max_localmem=lambda wildcards, resources: resources.mem_mb*resources.cpus_per_task/1000,
        samp_id=lambda wildcards: wildcards.pool, # Same as samp_id in rule "inp_cellranger_arc_count"
        outputdir=f"{config['cellranger_arc_count']['bams_dir']}{{pool}}",
        pipestance=f"{config['cellranger_arc_count']['pipestance_struct']}"

    output:
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['gex_bam']}",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['gex_bai']}",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['atac_bam']}",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['atac_bai']}",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['atac_fragments']}",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['filtered_h5_matrix']}",
        # expand(
        #     f"{config['cellranger_arc_count']['bams_dir']}/{{{pool}}}/filtered_feature_bc_matrix/{{name}}",
        #     name=["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"],
        # ),
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/filtered_feature_bc_matrix/features.tsv.gz",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/filtered_feature_bc_matrix/matrix.mtx.gz",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['pipestance_struct']}",
        outdir=protected(directory(
            f"{config['cellranger_arc_count']['bams_dir']}{{pool}}"
            ))

    log:
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/cellranger-arc-count.txt"
   
    resources:
        cpus_per_task=20, # For snakemake > v8
        mem_mb=allocate_mem_CAC,
        time_min=allocate_time_CAC,
        attempt=lambda wildcards, attempt: attempt

    # For snakemake < v8
    # threads: 20

    envmodules:
        f"{config['cellranger-arc_version']}"

    shell:
        """
        mkdir -p {params.outputdir}
        export MRO_DISK_SPACE_CHECK=disable
        loc_mem=$( sed r"s#\.0##g" <<< "{params.max_localmem}")
        cellranger-arc count --id={params.samp_id} --libraries={input[0]} \
        --reference={params.ref} --localcores={params.max_localcores} \
        --localmem=${{loc_mem}} &> {log}_{resources.attempt} && \
        rm -r {params.samp_id}/SC_ATAC_GEX_COUNTER_CS/ && \
        mv {params.samp_id}/outs/* {params.outputdir}/ && \
        mv {params.samp_id}/{params.pipestance} {params.outputdir}/
        """
