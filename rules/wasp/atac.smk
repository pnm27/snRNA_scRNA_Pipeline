# WASP of the ATAC bams (the ATAC part of multiome, scATAC): paired-end
# reads over SNPs are remapped with bwa-mem (as cellranger-arc), and the
# reads cellranger-arc marked as duplicates (0x400, per cell barcode) are
# removed after WASP ('drop_duplicates: marked'). Needs
# rules/wasp/common.smk. Settings: 'modules: wasp: atac'
# (lib/wasp_modality.py).

from lib.wasp_modality import bwa_index, bwa_extra, drop_duplicates


rule wasp_remap_atac:
    input:
        unpack(partial(wasp_input_bam, modality='ATAC')),
        unpack(partial(wasp_genotype_files, modality='ATAC')),
        checked = checked_inputs(config),

    output:
        **donor_paths(config, 'ATAC')

    params:
        **wasp_samples_param('ATAC'),
        bwa_index = lambda wc: bwa_index(config),
        bwa_extra = bwa_extra(config),
        drop_duplicates = lambda wc: "--drop_duplicates" if drop_duplicates(config) else "",
        threads = lambda wc, resources: resources.cpus_per_task,

    log:
        **donor_logs(config, 'ATAC')

    resources:
        **get_resources("wasp_remap_atac"),
        wasp_slots = 1,

    envmodules:
        *WASP_MODULES,
        "bwa",
        "samtools",

    shell:
        r"""
        set -euo pipefail
        work={params.work}
        bash helper_sh_scripts/wasp_donor.sh prepare --work "$work" --bam {input.bam} \
            --samples {params.samples:q} --snp_index {input.snp_index} \
            --snp_tab {input.snp_tab} --haps {input.haps} {params.paired} \
            --log {log.intersect}
        prefix="$work/$(basename {input.bam} .bam)"

        # Remap the pairs with the alleles swapped
        bwa mem -t {params.threads} {params.bwa_extra} {params.bwa_index} \
            "${{prefix}}.remap.fq1.gz" "${{prefix}}.remap.fq2.gz" 2> {log.remap} \
            | samtools view -@ {params.threads} -b -o "$work/remap.bam" -

        bash helper_sh_scripts/wasp_donor.sh finish --work "$work" --bam {input.bam} \
            --remapped "$work/remap.bam" --out {output.bam} \
            --threads {params.threads} --filter_log {log.filter} \
            --intersect_log {log.intersect} --flagstat_before {output.flagstat_before} \
            --flagstat_after {output.flagstat_after} \
            --flagstat_wasp {output.flagstat_wasp} {params.drop_duplicates} {params.paired} \
            --pool {wildcards.pool} --donor {params.donor} --final_donor {wildcards.final_donor} \
            --sample {params.sample_ids:q} --modality ATAC --stats {output.stats}
        """
