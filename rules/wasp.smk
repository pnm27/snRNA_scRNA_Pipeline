# WASP filtering of per-donor bams ('workflow: until: wasp', or 'start:
# wasp' for bams listed in the samples table), for sc/snRNA-seq and the
# RNA (GEX) part of multiome. See lib/wasp.py and docs/source/sub_workflows.md.
#
# Per donor bam (rule wasp_remap), as github.com/bmvdgeijn/WASP's mapping
# pipeline: find_intersecting_snps.py (reads over the donor's SNPs, the
# --samples IDs of the samples table), STAR remapping of those reads with
# the alleles swapped, filter_remapped_reads.py (keeps the reads that map
# back to the same place), merged with the reads over no SNP. Duplicates
# are kept (sc/snRNA-seq). WASP comes from its env module
# ('references: tools: wasp'), whose scripts are in $WASP_HOME.
#
# At most 'wasp_slots' jobs at a time (each holds one): the temporary bams
# are large, e.g. --resources wasp_slots=100.

from lib.wasp import (
    samples_table, starts_at_wasp, select_donors, missing_bams,
    donor_paths, donor_logs, work_dir, run_outputs, checked_inputs,
    genotype_files, genotype_sets, vcf_sets, vcf_files, converted_dir,
    remap_reference, wasp_config, H5_FILES,
)

WASP_TOOLS = config['references']['tools']
WASP_MODULES = [m for m in [WASP_TOOLS.get('wasp')] if m]
WASP_PAIRED = bool(wasp_config(config).get('paired_end'))


def wasp_row(wildcards) -> dict:
    row = samples_table(config).get((wildcards.pool, wildcards.donor))
    if row is None:
        raise ValueError(
            f"No row for pool {wildcards.pool}, donor {wildcards.donor} in "
            f"the WASP samples table {wasp_config(config)['samples']}")
    return row


def wasp_input_bam(wildcards) -> dict:
    """The donor's bam: from the samples table (start: wasp), else this
    run's split bam."""
    if starts_at_wasp(config):
        row = wasp_row(wildcards)
        return {'bam': row['bam'], 'bai': row['bai']}
    sb = config['modules']['split_bams']
    bam = f"{sb['split_bams_dir']}{sb['fold_struct_bams']}{wildcards.donor}.bam".format(
        pool=wildcards.pool, donor=wildcards.donor)
    return {'bam': bam, 'bai': f"{bam}.bai"}


def wasp_genotype_files(wildcards) -> dict:
    return genotype_files(config, wasp_row(wildcards)['genotype_set'])


_WASP_NOTED = set()


def wasp_donors() -> list:
    """(pool, donor) of every WASP job of the run: the rows of the samples
    table (start: wasp), else the split donors that are in the table."""
    table = samples_table(config)
    if starts_at_wasp(config):
        missing = missing_bams(table)
        if missing:
            raise ValueError(
                "'workflow: start: wasp' filters existing bams; not there:\n  "
                + "\n  ".join(missing))
        return sorted(table)
    pools = split_pools() if EXTERNAL else POOLS
    out, skipped = [], []
    for pool in pools:
        donors = get_real_donors({'pool': pool, 'modality': SPLIT_MODALITY})['donor']
        kept, missing = select_donors(table, pool, donors)
        out += [(pool, d) for d in kept]
        skipped += [f"{pool} {d}" for d in missing]
    if skipped and 'skipped' not in _WASP_NOTED:
        _WASP_NOTED.add('skipped')
        print("Warning: split donor(s) not in the WASP samples table "
              f"{wasp_config(config)['samples']}, not WASP-filtered: "
              + ", ".join(skipped), file=sys.stderr)
    if not out:
        raise ValueError(
            "No split donor is in the WASP samples table "
            f"{wasp_config(config)['samples']}: check its pool and donor columns")
    return out


def wasp_donor_stats(wildcards) -> list:
    return [donor_paths(config)['stats'].format(pool=p, donor=d) for p, d in wasp_donors()]


# Genotype sets given as VCF/BCF: WASP's h5 files, once per set ----------
rule wasp_snp2h5:
    input:
        vcf = lambda wc: vcf_files(config, wc.gset),
        table = wasp_config(config).get('samples') or [],
        chrom_sizes = lambda wc: wasp_config(config)['genotypes'][wc.gset]['chrom_sizes'],

    output:
        **{k: f"{converted_dir(config, '{gset}')}{f}" for k, f in H5_FILES.items()}

    # Only the sets given as VCF/BCF are made here
    wildcard_constraints:
        gset = "|".join(map(re.escape, vcf_sets(config))) or r"(?!)",

    params:
        sets = genotype_sets(config),
        need_bam = "--need_bam" if starts_at_wasp(config) else "",
        out_dir = lambda wc: converted_dir(config, wc.gset),

    log:
        f"{converted_dir(config, '{gset}')}snp2h5.log"

    resources:
        **get_resources("wasp_snp2h5")

    envmodules:
        *WASP_MODULES,
        "bcftools",

    shell:
        r"""
        set -euo pipefail
        : "${{WASP_HOME:?WASP_HOME isn't set: load the WASP module (references: tools: wasp)}}"
        out={params.out_dir}
        mkdir -p "$out/vcf"
        python helper_py_scripts/wasp_check_inputs.py set-samples --table {input.table} \
            --set {wildcards.gset} --sets {params.sets} {params.need_bam} > "$out/samples.txt"
        # One file per chromosome (of chrom_sizes), with the set's samples
        # and the sites variable among them (the vcfs must be indexed)
        i=0
        for vcf in {input.vcf}; do
            i=$((i + 1))
            for chrom in $(bcftools index -s "$vcf" | cut -f1); do
                awk -v c="$chrom" '$1 == c {{f=1}} END {{exit !f}}' {input.chrom_sizes} || continue
                bcftools view -S "$out/samples.txt" -q 0.000001:minor -r "$chrom" "$vcf" \
                    -Oz -o "$out/vcf/${{chrom}}_${{i}}.vcf.gz"
            done
        done &> {log}
        if ! ls "$out"/vcf/*.vcf.gz &> /dev/null; then
            echo "No chromosome of {input.chrom_sizes} in {input.vcf}" | tee -a {log} >&2
            exit 1
        fi
        $WASP_HOME/snp2h5/snp2h5 --chrom {input.chrom_sizes} --format vcf \
            --snp_index {output.snp_index} --snp_tab {output.snp_tab} \
            --haplotype {output.haps} "$out"/vcf/*.vcf.gz &>> {log}
        rm -rf "$out/vcf"
        """


# The samples table's IDs are in their sets' haps.h5, once per run --------
rule wasp_check_inputs:
    input:
        table = wasp_config(config).get('samples') or [],
        h5 = [f for s in genotype_sets(config) for f in genotype_files(config, s).values()],

    output:
        checked_inputs(config)

    params:
        sets = " ".join(
            f"--set {s}={genotype_files(config, s)['haps']},{genotype_files(config, s)['snp_tab']}"
            for s in genotype_sets(config)),
        need_bam = "--need_bam" if starts_at_wasp(config) else "",

    resources:
        **get_resources("wasp_check_inputs")

    envmodules:
        *WASP_MODULES,

    shell:
        r"""
        python helper_py_scripts/wasp_check_inputs.py samples --table {input.table} \
            {params.sets} {params.need_bam} -o {output}
        """


# One donor bam ------------------------------------------------------------
rule wasp_remap:
    input:
        unpack(wasp_input_bam),
        unpack(wasp_genotype_files),
        checked = checked_inputs(config),

    output:
        **donor_paths(config)

    params:
        samples = lambda wc: wasp_row(wc)['samples'],
        sample_ids = lambda wc: ";".join(wasp_row(wc)['samples']),
        work = lambda wc: work_dir(config).format(pool=wc.pool, donor=wc.donor),
        genome = lambda wc: remap_reference(config)[0],
        gtf = lambda wc: (
            f"--sjdbGTFfile {remap_reference(config)[1]}" if remap_reference(config)[1] else ""),
        star_extra = (wasp_config(config).get('remap') or {}).get('extra_params') or "",
        paired = "--is_paired_end" if WASP_PAIRED else "",
        paired_stats = "--paired_end" if WASP_PAIRED else "",
        # (shell variables, so a function: no wildcard expansion)
        fastqs = lambda wc: (
            "${prefix}.remap.fq1.gz ${prefix}.remap.fq2.gz" if WASP_PAIRED
            else "${prefix}.remap.fq.gz"),
        threads = lambda wc, resources: resources.cpus_per_task,

    log:
        **donor_logs(config)

    resources:
        **get_resources("wasp_remap"),
        wasp_slots = 1,

    envmodules:
        *WASP_MODULES,
        config['references']['tools']['star'],
        "samtools",

    shell:
        r"""
        set -euo pipefail
        : "${{WASP_HOME:?WASP_HOME isn't set: load the WASP module (references: tools: wasp)}}"
        work={params.work}
        rm -rf "$work"
        mkdir -p "$work"
        printf '%s\n' {params.samples:q} > "$work/samples.txt"

        # The bam and the genotypes name chromosomes the same way
        samtools view -H {input.bam} \
            | awk -F'\t' '$1 == "@SQ" {{for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) print substr($i, 4)}}' \
            > "$work/contigs.txt"
        python helper_py_scripts/wasp_check_inputs.py contigs --contigs "$work/contigs.txt" \
            --snp_tab {input.snp_tab} --bam {input.bam}

        # Reads over the donor's SNPs (to remap) and the others (kept)
        python $WASP_HOME/mapping/find_intersecting_snps.py {params.paired} --is_sorted \
            --snp_index {input.snp_index} --snp_tab {input.snp_tab} \
            --haplotype {input.haps} --samples "$work/samples.txt" \
            --output_dir "$work" {input.bam} &> {log.intersect}
        prefix="$work/$(basename {input.bam} .bam)"

        # Remap them with the alleles swapped
        STAR --runThreadN {params.threads} --genomeDir {params.genome} {params.gtf} \
            --readFilesIn {params.fastqs} --readFilesCommand zcat \
            --twopassMode Basic --limitOutSJcollapsed 2000000 --limitSjdbInsertNsj 2000000 \
            --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD \
            --limitBAMsortRAM 0 --outBAMsortingBinsN 50 \
            --outFileNamePrefix "$work/remap." {params.star_extra} &> {log.star}
        cat "$work/remap.Log.final.out" >> {log.star}

        # Keep the reads that map back to the same place, with those over no SNP
        python $WASP_HOME/mapping/filter_remapped_reads.py "${{prefix}}.to.remap.bam" \
            "$work/remap.Aligned.sortedByCoord.out.bam" "$work/remap.keep.bam" &> {log.filter}
        samtools merge -f -@ {params.threads} "$work/merged.bam" \
            "${{prefix}}.keep.bam" "$work/remap.keep.bam"
        samtools sort -@ {params.threads} -T "$work/sort" -o {output.bam} "$work/merged.bam"
        samtools index -@ {params.threads} {output.bam}

        samtools flagstat -@ {params.threads} {input.bam} > {output.flagstat_before}
        samtools flagstat -@ {params.threads} {output.bam} > {output.flagstat_after}
        python helper_py_scripts/wasp_stats.py parse --pool {wildcards.pool} \
            --donor {wildcards.donor} --sample {params.sample_ids:q} {params.paired_stats} \
            --intersect_log {log.intersect} --filter_log {log.filter} \
            --flagstat_before {output.flagstat_before} --flagstat_after {output.flagstat_after} \
            -o {output.stats}
        rm -rf "$work"
        """


# The run's stats: one row per donor bam, a summary and plots -------------
rule wasp_stats:
    input:
        wasp_donor_stats

    output:
        **run_outputs(config)

    resources:
        **get_resources("wasp_stats")

    conda: "../envs/basic_sctools.yaml"

    shell:
        r"""
        python3 helper_py_scripts/wasp_stats.py summarize {input} \
            --stats {output.stats} --summary {output.summary} --plots {output.plots}
        """
