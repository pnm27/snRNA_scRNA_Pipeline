# WASP of the cDNA bams (sc/snRNA-seq, the RNA part of multiome): the
# reads over SNPs are remapped with STAR (module 'references: tools:
# wasp_star', else 'star'; a gzipped gtf is unzipped for it, then
# removed); duplicates are kept. Its own job, between the donor's
# wasp_prepare_rna and wasp_finish_rna (rules/wasp/common.smk), so a failed
# STAR is rerun alone. Settings: 'modules: wasp: rna' (lib/wasp_modality.py).

from lib.wasp_modality import remap_reference, rna_config, remap_star_module


rule wasp_remap_rna:
    input:
        fastqs = remap_fastqs(config, 'cDNA', paired_end(config, 'cDNA')),

    output:
        bam = temp(remapped_bam(config, 'cDNA')),

    params:
        scratch = remap_scratch(config, 'cDNA'),
        genome = lambda wc: remap_reference(config)[0],
        # a gzipped gtf is unzipped into the scratch dir for STAR, and removed after
        gtf = lambda wc: remap_reference(config)[1] or "",
        star_extra = (rna_config(config).get('remap') or {}).get('extra_params') or "",
        threads = lambda wc, resources: resources.cpus_per_task*3,

    log:
        donor_logs(config, 'cDNA')['remap']

    resources:
        **get_resources("wasp_remap_rna"),
        wasp_slots = 1,

    envmodules:
        remap_star_module(config),  # references: tools: wasp_star, else star

    shell:
        r"""
        set -euo pipefail
        # STAR's own files; a new start: nothing left of an earlier try
        scratch={params.scratch}
        scratch=${{scratch%/}}
        rm -rf "$scratch"
        mkdir -p "$scratch"

        # The gtf (STAR reads it plain): a gzipped one unzipped into the
        # scratch dir, removed when STAR is done, or the job fails
        gtf={params.gtf:q}
        gtf_arg=""
        if [ -n "$gtf" ]; then
            case "$gtf" in
                *.gz)
                    trap 'rm -f "$scratch/remap.gtf"' EXIT
                    gzip -cd "$gtf" > "$scratch/remap.gtf"
                    gtf="$scratch/remap.gtf" ;;
            esac
            gtf_arg="--sjdbGTFfile $gtf"
        fi

        # Remap the reads over SNPs, with the alleles swapped
        STAR --runThreadN {params.threads} --genomeDir {params.genome} $gtf_arg \
            --readFilesIn {input.fastqs} --readFilesCommand zcat \
            --twopassMode Basic --limitOutSJcollapsed 2000000 --limitSjdbInsertNsj 2000000 \
            --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD \
            --limitBAMsortRAM 0 --outBAMsortingBinsN 50 \
            --outFileNamePrefix "$scratch/" {params.star_extra} &> {log}
        cat "$scratch/Log.final.out" >> {log}
        mv "$scratch/Aligned.sortedByCoord.out.bam" {output.bam}
        rm -rf "$scratch"
        """
