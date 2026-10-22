from functools import partial
from lib.io import get_bam_inputs


rule qtltools_mbv:
    input:
        bam=partial(get_bam_inputs, config=config),
        ref_snps=config['modules']['identify_swaps']['ref_vcf']

    output:
        (
            f"{config['modules']['identify_swaps']['mbv_out_dir']}"
            f"{config['modules']['identify_swaps']['fold_struct']}"
            f"{config['modules']['identify_swaps']['mbv_suffix']}"
        )
        

    params:
        extra_params=config['modules']['identify_swaps']['mbv_extra_opt']

    resources:
        **get_resources("qtltools_mbv")

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