# WASP of the cDNA bams (sc/snRNA-seq, the RNA part of multiome): the
# reads over SNPs are remapped with STAR; duplicates are kept. Needs
# rules/wasp/common.smk. Settings: 'modules: wasp: rna'
# (lib/wasp_modality.py).

from lib.wasp_modality import remap_reference, rna_config


rule wasp_remap_rna:
    input:
        unpack(partial(wasp_input_bam, modality='cDNA')),
        unpack(partial(wasp_genotype_files, modality='cDNA')),
        checked = checked_inputs(config),

    output:
        **donor_paths(config, 'cDNA')

    params:
        **wasp_samples_param('cDNA'),
        genome = lambda wc: remap_reference(config)[0],
        gtf = lambda wc: (
            f"--sjdbGTFfile {remap_reference(config)[1]}" if remap_reference(config)[1] else ""),
        star_extra = (rna_config(config).get('remap') or {}).get('extra_params') or "",
        # (shell variables, so a function: no wildcard expansion)
        fastqs = lambda wc: (
            "${prefix}.remap.fq1.gz ${prefix}.remap.fq2.gz" if paired_end(config, 'cDNA')
            else "${prefix}.remap.fq.gz"),
        threads = lambda wc, resources: resources.cpus_per_task,

    log:
        **donor_logs(config, 'cDNA')

    resources:
        **get_resources("wasp_remap_rna"),
        wasp_slots = 1,

    envmodules:
        *WASP_MODULES,
        config['references']['tools']['star'],
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

        # Remap them with the alleles swapped
        STAR --runThreadN {params.threads} --genomeDir {params.genome} {params.gtf} \
            --readFilesIn {params.fastqs} --readFilesCommand zcat \
            --twopassMode Basic --limitOutSJcollapsed 2000000 --limitSjdbInsertNsj 2000000 \
            --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD \
            --limitBAMsortRAM 0 --outBAMsortingBinsN 50 \
            --outFileNamePrefix "$work/remap." {params.star_extra} &> {log.remap}
        cat "$work/remap.Log.final.out" >> {log.remap}

        bash helper_sh_scripts/wasp_donor.sh finish --work "$work" --bam {input.bam} \
            --remapped "$work/remap.Aligned.sortedByCoord.out.bam" --out {output.bam} \
            --threads {params.threads} --filter_log {log.filter} \
            --intersect_log {log.intersect} --flagstat_before {output.flagstat_before} \
            --flagstat_after {output.flagstat_after} {params.paired} \
            --pool {wildcards.pool} --donor {wildcards.donor} \
            --sample {params.sample_ids:q} --modality cDNA --stats {output.stats}
        """
