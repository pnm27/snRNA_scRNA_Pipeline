import os, glob2, re
from functools import partial
from lib.io import get_fastqs, get_r1_fastqs, get_r2_fastqs
from lib.aligner import bam_paths, out_prefix


rule inp_cellranger_arc_count:
    input:
        partial(get_fastqs, config=config)

    params:
        samp_id=lambda wildcards: wildcards.pool # Same as samp_id in rule "cellranger_arc_count"

    output:
        "{pool}_metadata.csv"

    resources:
        **get_resources("inp_cellranger_arc_count")

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
        partial(get_fastqs, config=config)

    params:
        ref=config['cellranger_arc_count']['ref']['genome'],
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
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/filtered_feature_bc_matrix/features.tsv.gz",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/filtered_feature_bc_matrix/matrix.mtx.gz",
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{config['cellranger_arc_count']['pipestance_struct']}",
        outdir=directory(
            f"{config['cellranger_arc_count']['bams_dir']}{{pool}}"
            )

    log:
        f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/cellranger-arc-count.txt"
    
    # Restrict wildcards so that ATAC and cDNA filt or subset bams 
    # (from split_bams.smk) won't get triggerred
    # wildcard_constraints:
    #     pool=r"[^/]" #WILDCARDS

    resources:
        **get_resources("cellranger_arc_count"),
        attempt=lambda wildcards, attempt: attempt

    envmodules:
        f"{config['cellranger-arc_version']}"

    shell:
        r"""
        mkdir -p {params.outputdir}
        export MRO_DISK_SPACE_CHECK=disable
        loc_mem=$( sed "s#\.0##g" <<< "{params.max_localmem}")
        cellranger-arc count --id={params.samp_id} --libraries={input[0]} \
        --reference={params.ref} --localcores={params.max_localcores} \
        --localmem=${{loc_mem}} &> {log}_{resources.attempt} && \
        rm -rf {params.samp_id}/SC_ATAC_GEX_COUNTER_CS/ && \
        rm -rf {params.samp_id}/SC_RNA_COUNTER_CS/ && \
        mv {params.samp_id}/outs/* {params.outputdir}/ && \
        mv {params.samp_id}/{params.pipestance} {params.outputdir}/
        """


# cellranger count (GEX-only runs, last_step 'cellranger_*') ----------------

def _fastq_sample(fastq: str) -> str:
    """cellranger's --sample of a fastq (its name before _S<n>_L<nnn>_R1)."""
    name = os.path.basename(fastq)
    m = re.match(r"(.+?)_S\d+_L\d{3}_R[12]_\d{3}\.fastq(\.gz)?$", name)
    if m is None:
        raise ValueError(
            f"Can't get cellranger's sample name from '{name}' (expected "
            "<sample>_S<n>_L<nnn>_R1_001.fastq.gz)"
        )
    return m.group(1)


rule cellranger_count:
    """
    Outputs of 'cellranger count' are moved from <pool>/outs/ to
    'cellranger_count: bams_dir'<pool>/.
    """
    input:
        R1=partial(get_r1_fastqs, config=config),
        R2=partial(get_r2_fastqs, config=config)

    params:
        ref=config['cellranger_count']['ref']['genome'],
        fastqs_dir=lambda wildcards, input: os.path.dirname(input.R1[0]),
        sample=lambda wildcards, input: ",".join(
            sorted({_fastq_sample(f) for f in input.R1})
        ),
        include_introns=str(config['cellranger_count']['include_introns']).lower(),
        extra=config['cellranger_count']['extra_params'] or "",
        max_localcores=lambda wildcards, resources: resources.cpus_per_task,
        max_localmem=lambda wildcards, resources: resources.mem_mb*resources.cpus_per_task//1000,
        outputdir=out_prefix(config, 'cellranger')

    output:
        **bam_paths(config, 'cDNA', aligner='cellranger'),
        matrix=f"{out_prefix(config, 'cellranger')}filtered_feature_bc_matrix/matrix.mtx.gz",
        features=f"{out_prefix(config, 'cellranger')}filtered_feature_bc_matrix/features.tsv.gz",
        barcodes=f"{out_prefix(config, 'cellranger')}filtered_feature_bc_matrix/barcodes.tsv.gz",
        metrics=f"{out_prefix(config, 'cellranger')}metrics_summary.csv"

    log:
        f"{out_prefix(config, 'cellranger')}cellranger-count.txt"

    resources:
        **get_resources("cellranger_count"),
        attempt=lambda wildcards, attempt: attempt

    envmodules:
        f"{config['cellranger_version']}"

    shell:
        r"""
        mkdir -p {params.outputdir}
        rm -rf {wildcards.pool}
        cellranger count --id={wildcards.pool} --transcriptome={params.ref}             --fastqs={params.fastqs_dir} --sample={params.sample}             --create-bam=true --include-introns={params.include_introns}             --localcores={params.max_localcores} --localmem={params.max_localmem}             {params.extra} &> {log}_{resources.attempt} &&         mv {wildcards.pool}/outs/* {params.outputdir}/ &&         rm -rf {wildcards.pool}
        """
