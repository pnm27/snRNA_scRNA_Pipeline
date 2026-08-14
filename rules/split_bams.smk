# import os, uuid
import pandas as pd
from lib.io import (
    get_inp_splitBam, 
    get_bam,
    get_bam_to_split
)
from lib.params import (
    get_params_create_inp_splitBams,
    get_mito, 
    subset_to_chr, 
    get_mito_file,
)
from functools import partial


# Resource Allocation ------------------
def allocate_mem_FCB(wildcards, attempt):
    return 1000+200*(attempt-1)


def allocate_time_FCB(wildcards, attempt):
    return 20+20*(attempt-1)


def allocate_mem_SB(wildcards, attempt):
    return 250+100*(attempt-1)

# --------------------------------------

# FIX THE INPUT: INPUT SHOULD TAKE IN EITHER VIREO OR CALICO_SOLO'S H5AD AS INPUT, AS REQUIRED.
checkpoint create_inp_splitBams:
    input:
        partial(get_inp_splitBam, config=config)

    output:
        (
            f"{config['split_bams_pipeline']['inp_split_bams_dir']}"
            f"{config['fold_struct_bam_split1']}_bc_hash.txt"
        )

    params:
        # DEPRACATED
        # conv=config['split_bams_pipeline']['donor_name_converter']['file'],
        # from_col=config['split_bams_pipeline']['donor_name_converter']['from_column'],
        # to_col=config['split_bams_pipeline']['donor_name_converter']['to_column'],
        # demux_method=config['split_bams_pipeline']['split_by']['demux'],
        # inp_ext=config['split_bams_pipeline']['split_by']['input'],
        # h5ad_col=config['split_bams_pipeline']['split_by']['column'],
        extra = partial(get_params_create_inp_splitBams, config=config),
        script="helper_py_scripts/create_inp_splitBam_consolidated.py"

    # For snakemake < v8
    # threads: 2

    resources:
        cpus_per_task=2, # For snakemake > v8
        mem_mb=1000,
        time_min=10

    shell:
        """
        python3 {params.script} {input} {output} {params.extra}

        sleep 30
        """

# Generalize Output
rule bamfilt_by_CB:
    input:
        unpack(partial(get_bam, config=config)),
        hash_file=(
            f"{config['split_bams_pipeline']['inp_split_bams_dir']}"
            f"{config['fold_struct_bam_split1']}_bc_hash.txt"
            )

    output:
        bam=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['filt_bam']}"
        ),
        bai=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['filt_bam']}.bai"
        )

    params:
        temp_bc=lambda wc: f"{config['split_bams_pipeline']['sort_temp_dir']}{{pool}}_bc.txt",
        threads=lambda _, resources: max(resources.cpus_per_task*6, 10)


    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=250,
        time_min=300

    conda: "../envs/pysam.yaml"

    envmodules:
        "samtools"

    shell:
        r"""
        mkdir -p {config[split_bams_pipeline][sort_temp_dir]}
        cut -f2 <(tail -n +2 {input.hash_file}) > {params.temp_bc}
        samtools view -@ {params.threads} -q 255 -D CB:{params.temp_bc} \
            {input.bam} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        rm {params.temp_bc}
        sleep 10
        """


rule filt_chr_bams:
    input:
        unpack(partial(get_bam, config=config))

    output:
        bam=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['short_bam']}"
        ), # generalize this
        bai=(
            f"{config['STARsolo_pipeline']['bams_dir']}"
            f"{config['fold_struct']}"
            f"{config['split_bams_pipeline']['short_bam']}.bai"
        )

    params:
        sub_chr=partial(subset_to_chr, config=config),
        threads=lambda _, resources: max(resources.cpus_per_task*6, 10)

    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_FCB,
        time_min=allocate_time_FCB

    conda: "../envs/pysam.yaml"

    envmodules:
        "samtools"
        
    shell:
        """
        samtools view -@ {params.threads} {input.bam} {params.sub_chr} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        """


def _bam_pattern(config, modality_key: str) -> str:
    """Return a wildcard-pattern path for a given modality bam key."""
    bam_file = config['cellranger_arc_count'][modality_key]
    short    = config['split_bams_pipeline']['short_bam']
    base     = (
        f"{config['cellranger_arc_count']['bams_dir']}"
        f"{{pool}}/{bam_file}"
    ).replace('.bam', short)
    return base


def _output_fcb_multiome(config) -> tuple[str, str]:
    """
    Return (bam_pattern, bai_pattern) for multiome, driven by {modality} wildcard.
    The {modality} wildcard in the pattern lets Snakemake resolve it at runtime.
    """
    gex_bam  = _bam_pattern(config, 'gex_bam')
    atac_bam = _bam_pattern(config, 'atac_bam')

    # Use a Snakemake wildcard-keyed lookup at rule-match time via a
    # dispatch pattern encoded directly into the path using {modality}
    # We can't branch on modality here, so we rely on the {modality}
    # wildcard constraint on the rule to route correctly — two rules,
    # one per modality, or a shared pattern if paths can encode modality.
    #
    # Since gex and atac have different filenames, we need two rules.
    # Return both as a dict for callers to pick from.
    return {
        "cdna": (gex_bam,  gex_bam  + '.bai'),
        "atac": (atac_bam, atac_bam + '.bai'),
    }


_out_fcb_multiome = _output_fcb_multiome(config)

rule filt_chr_bams_multiome_cdna:
    wildcard_constraints:
        modality="cDNA"
    input:
        # unpack(partial(get_bam, config=config))
        unpack(lambda wildcards: get_bam(
            SimpleNamespace(**vars(wildcards), modality="cDNA"), config
        ))
    output:
        bam=_out_fcb_multiome["cdna"][0],
        bai=_out_fcb_multiome["cdna"][1]
    params:
        sub_chr=partial(subset_to_chr, config=config),
        threads=lambda _, resources: max(resources.cpus_per_task * 6, 10)
    resources:
        cpus_per_task=1,
        mem_mb=allocate_mem_FCB,
        time_min=allocate_time_FCB
    conda: "../envs/pysam.yaml"
    envmodules:
        "samtools"
    shell:
        """
        samtools view -@ {params.threads} {input.bam} {params.sub_chr} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        """


rule filt_chr_bams_multiome_atac:
    wildcard_constraints:
        modality="ATAC"
    input:
        # unpack(partial(get_bam, config=config))
        unpack(lambda wildcards: get_bam(
            SimpleNamespace(**vars(wildcards), modality="ATAC"), config
        ))
    output:
        bam=_out_fcb_multiome["atac"][0],
        bai=_out_fcb_multiome["atac"][1]
    params:
        sub_chr=partial(subset_to_chr, config=config),
        threads=lambda _, resources: max(resources.cpus_per_task * 6, 10)
    resources:
        cpus_per_task=1,
        mem_mb=allocate_mem_FCB,
        time_min=allocate_time_FCB
    conda: "../envs/pysam.yaml"
    envmodules:
        "samtools"
    shell:
        """
        samtools view -@ {params.threads} {input.bam} {params.sub_chr} -bho {output.bam}
        samtools index -@ {params.threads} {output.bam}
        """


rule create_bed:
    input:
        config['STARsolo_pipeline']['genome_pick']['fasta']

    params:
        chr_prefix= lambda _: (
            '' if config['chr_prefix'] is None 
            else config['chr_prefix']
        )

    output:
        config['reg_chr_bed']

    shell:
        r"""
        if [ ! -f "{output}" ] ; then 
            grep -E "^>{params.chr_prefix}[0-9]+|X|Y|MT" "{input}" | \
            awk 'BEGIN{{OFS="\t"}}{{split($3,a,":");gsub(">", "", $1);printf("%s\\t0\\t%s\\n",\$1,a[5]);}}' > "{output}"
        fi
        sleep 10
        """


def get_real_donors(wildcards):
    import pandas as pd
    from itertools import repeat

    ckpt = checkpoints.create_inp_splitBams.get(**wildcards)
    demux_out = ckpt.output[0]
    min_cells_threshold = config['split_bams_pipeline']['min_cells_per_donor']

    df = pd.read_csv(demux_out, sep="\t")

    # filter donors with meaningful data
    df = df.groupby('Subj_ID')['barcodes'].count()
    donors = df[df>min_cells_threshold].index.tolist()

    return {'pool': list(repeat(wildcards['pool'], len(donors))) , 'donor': donors}

# ---------------------------------------------------------------------------
# Rule: one job per donor, both gt_check modes
# ---------------------------------------------------------------------------

rule split_bams:
    input:
        unpack(partial(get_bam_to_split, config=config)),
        barcodes_vs_donor = (
            f"{config['split_bams_pipeline']['inp_split_bams_dir']}"
            f"{config['fold_struct_bam_split1']}_bc_hash.txt"
        )
    output:
        bam = (
            f"{config['split_bams_pipeline']['split_bams_dir']}"
            f"{config['fold_struct_bam_split2']}"
            "{donor}.bam"
        ),
        bai = (
            f"{config['split_bams_pipeline']['split_bams_dir']}"
            f"{config['fold_struct_bam_split2']}"
            "{donor}.bam.bai"
        )
    params:
        temp_dir      = (
            f"{config['split_bams_pipeline']['new_temp_dir']}"
            f"{config['fold_struct_bam_split2']}"
        ),
        gt_check      = config["gt_check"],
        chr_mito      = partial(get_mito, config=config),
        mito_info     = partial(get_mito_file, config=config),
        filter_bed    = lambda _: (
            None if config["gt_check"]
            else config.get("reg_chr_bed")
        ),
        script        = "helper_sh_scripts/create_per_donor_bams.bash"
    resources:
        cpus_per_task = 1,
        mem_mb        = allocate_mem_SB,
        time_min      = 30
    log:
        (
            f"{config["split_bams_pipeline"]["per_donor_split_log_dir"]}"
            f"{config['fold_struct_bam_split2']}"
            "{donor}.log",
        )
    shell:
        """
        # Build the optional-argument string based on mode.
        # gt_check mode: no BED or mito refinement needed.
        # full-split mode: pass BED and/or mito file if configured.

        EXTRA_ARGS=""

        if [[ "{params.gt_check}" != "True" ]]; then
            if [[ "{params.filter_bed}" != "None" && "{params.filter_bed}" != "False" ]]; then
                EXTRA_ARGS+=" --bed {params.filter_bed}"
            fi
            if [[ "{params.mito_info}" != "None" && "{params.mito_info}" != "False" ]]; then
                EXTRA_ARGS+=" --mito_file {params.mito_info} --mito_pref {params.chr_mito}"
            fi
        fi

        bash {params.script} \
            --donor        {wildcards.donor}         \
            --hash_file    {input.barcodes_vs_donor} \
            --temp_dir     {params.temp_dir}         \
            --out_dir      $(dirname {output.bam})   \
            --bam          {input.bam}          \
            ${{EXTRA_ARGS}}                          \
            > {log} 2>&1
        """
