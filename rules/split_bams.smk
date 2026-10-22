# import os, uuid
import pandas as pd
from lib.io import (
    split_tool,
    get_bam,
    bam_to_split,
)
from lib.config_models import demux_run, final_paths, raw_output
from lib.params import (
    get_params_create_inp_splitBams,
    get_mito, 
    subset_to_chr, 
    get_mito_file,
)
from functools import partial
from lib.modality import bam_modalities, is_multiome
from lib.aligner import bam_paths


# Donor classification the bams are split by ('split_by: input'): the
# output of the demultiplexing tool itself ('raw'; multiome: vireo's cDNA
# output, as only cDNA bams are split for now), its h5ad in this run
# ('h5ad'), or the final h5ad of the pool, e.g. after the swap correction
checkpoint create_inp_splitBams:
    input:
        branch(
            lambda _: config['split_bams_pipeline']['split_by']['input'].lower(),
            cases={
                'raw': lambda _: raw_output(config, split_tool(config)).replace('{modality}', 'cDNA'),
                'h5ad': lambda _: demux_run(config, split_tool(config), config['run_tag']).h5ad,
                'final': lambda _: final_paths(config).h5ad,
            },
        )

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

    resources:
        **get_resources("create_inp_splitBams")

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

    # cDNA only for now (see TODO: split bams for multiome)
    output:
        **bam_paths(config, 'cDNA', config['split_bams_pipeline']['filt_bam'])

    params:
        temp_bc=lambda wc: f"{config['split_bams_pipeline']['sort_temp_dir']}{wc.pool}_bc.txt",
        threads=lambda _, resources: max(resources.cpus_per_task*6, 10)


    resources:
        **get_resources("bamfilt_by_CB")

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


# One rule per modality with a bam ('filt_chr_bams' for STARsolo runs)
for _modality in bam_modalities(config):
    _name = (
        f"filt_chr_bams_multiome_{_modality.lower()}" if is_multiome(config)
        else "filt_chr_bams"
    )

    rule:
        name: _name

        input:
            **bam_paths(config, _modality)

        output:
            **bam_paths(config, _modality, config['split_bams_pipeline']['short_bam'])

        params:
            sub_chr=partial(subset_to_chr, config=config),
            threads=lambda _, resources: max(resources.cpus_per_task*6, 10)

        resources:
            **get_resources(_name)

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
        unpack(branch(
            partial(bam_to_split, config=config),
            cases={
                'whole': partial(get_bam, config=config),
                'chr_subset': partial(get_bam, config=config,
                    suffix=config['split_bams_pipeline']['short_bam']),
                'cell_filtered': partial(get_bam, config=config,
                    suffix=config['split_bams_pipeline']['filt_bam']),
            },
        )),
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
        **get_resources("split_bams")
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
