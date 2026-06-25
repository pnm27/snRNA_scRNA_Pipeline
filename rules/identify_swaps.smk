from functools import partial
from lib.io import get_bam_inputs


# Resource Allocation ------------------
def allocate_mem_QM(wildcards):
    return 5000


def allocate_time_QM(wildcards, attempt):
    return 60 + 20*(attempt-1)

# --------------------------------------

rule qtltools_mbv:
    input:
        bam=partial(get_bam_inputs, config=config),
        ref_snps=config['identify_swaps']['ref_vcf']

    output:
        (
            f"{config['identify_swaps']['mbv_out_dir']}"
            f"{config['fold_struct_swaps_check']}"
            f"{config['identify_swaps']['mbv_suffix']}"
        )
        

    params:
        extra_params=config['identify_swaps']['mbv_extra_opt']

    # For snakemake < v8
    # threads: 1

    resources:
        cpus_per_task=1, # For snakemake > v8
        mem_mb=allocate_mem_QM,
        time_min=allocate_time_QM

    envmodules:
        "qtltools/1.3"

    # Use one of the below
    wildcard_constraints:
        # donor=r"(?<=cDNA_|ATAC_).+" #WILDCARDS # For multiome
        # donor=r"(?:.*)(?<=\/)([^/_]*)", #WILDCARDS # Match everything except the last '/', donor var should not contain '/' and '_'
        # donor=r"[^.]+",
        pool=r"[^_]+"

    shell:
        """
        if [[ "{params.extra_params}" == "None" ]]; then
            QTLtools mbv --bam {input.bam} --out {output} --vcf {input.ref_snps}
        else
            QTLtools mbv --bam {input.bam} --out {output} --vcf {input.ref_snps} {params.extra_params}
        fi
        
        """