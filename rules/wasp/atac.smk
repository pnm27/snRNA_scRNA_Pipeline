# WASP of the ATAC bams (the ATAC part of multiome, scATAC): paired-end
# reads over SNPs are remapped with bwa-mem (as cellranger-arc), and the
# reads cellranger-arc marked as duplicates (0x400, per cell barcode) are
# removed after WASP ('drop_duplicates: marked', in wasp_finish_atac). Its
# own job, between the donor's wasp_prepare_atac and wasp_finish_atac
# (rules/wasp/common.smk), so a failed bwa is rerun alone. Settings:
# 'modules: wasp: atac' (lib/wasp_modality.py).

from lib.wasp_modality import bwa_index, bwa_extra


rule wasp_remap_atac:
    input:
        fastqs = remap_fastqs(config, 'ATAC', paired_end(config, 'ATAC')),

    output:
        bam = temp(remapped_bam(config, 'ATAC')),

    params:
        bwa_index = lambda wc: bwa_index(config),
        bwa_extra = bwa_extra(config),
        threads = lambda wc, resources: resources.cpus_per_task*3,

    log:
        donor_logs(config, 'ATAC')['remap']

    resources:
        **get_resources("wasp_remap_atac"),
        wasp_slots = 1,

    envmodules:
        "bwa",
        "samtools",

    shell:
        r"""
        set -euo pipefail
        # Remap the pairs over SNPs, with the alleles swapped
        bwa mem -t {params.threads} {params.bwa_extra} {params.bwa_index} \
            {input.fastqs} 2> {log} \
            | samtools view -@ {params.threads} -b -o {output.bam} -
        """
