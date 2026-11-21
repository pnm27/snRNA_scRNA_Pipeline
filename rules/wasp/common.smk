# WASP filtering of per-donor bams ('workflow: until: wasp', or 'start:
# wasp' for bams listed in the samples table): what every modality shares.
# See lib/wasp.py and docs/source/sub_workflows.md.
#
#   common.smk  genotypes (snp2h5), input checks, the run's stats, and the
#               input functions of the per-donor jobs
#   rna.smk     cDNA (sc/snRNA-seq, multiome GEX): STAR remapping
#   atac.smk    ATAC (multiome ATAC, scATAC): bwa-mem remapping, duplicates
#               marked by cellranger-arc removed
#
# A pipeline of one modality includes common.smk and its own file. Per donor
# bam, as github.com/bmvdgeijn/WASP's mapping pipeline:
# find_intersecting_snps.py (reads over the donor's SNPs, the --samples IDs
# of the samples table), remapping of those reads with the alleles swapped,
# filter_remapped_reads.py (keeps the reads that map back to the same
# place), merged with the reads over no SNP (helper_sh_scripts/wasp_donor.sh
# runs all but the remapping). WASP comes from its env module
# ('references: tools: wasp'), whose scripts are in $WASP_HOME.
#
# At most 'wasp_slots' jobs at a time (each holds one): the temporary bams
# are large, e.g. --resources wasp_slots=100.

from lib.wasp import (
    samples_table, starts_at_wasp, select_donors, missing_bams, table_row, donor_of,
    donor_paths, donor_logs, work_dir, run_outputs, modality_dir, checked_inputs,
    genotype_files, genotype_sets, vcf_sets, vcf_files, converted_dir,
    wasp_config, sample_names, H5_FILES,
)
from lib.wasp_modality import paired_end
from lib.utils import split_bams_fold

ALL_TOOLS = config['references']['tools']
WASP_MODULES = [m for m in [ALL_TOOLS.get('wasp')] if m]
WASP_NAMES = sample_names(config)
WASP_NAMES_ARG = f"--names {WASP_NAMES}" if WASP_NAMES else ""


def wasp_donor(wildcards, modality: str) -> str:
    """The donor (as split; the samples table's donor) of a job's final donor."""
    donor = donor_of(samples_table(config), wildcards.pool, wildcards.final_donor, modality)
    if donor is None:
        raise ValueError(
            f"No row for pool {wildcards.pool}, final donor {wildcards.final_donor} "
            f"({modality}) in the WASP samples table {wasp_config(config)['samples']}")
    return donor


def wasp_row(wildcards, modality: str) -> dict:
    return table_row(samples_table(config), wildcards.pool,
                     wasp_donor(wildcards, modality), modality)


def wasp_input_bam(wildcards, modality: str) -> dict:
    """The donor's bam of 'modality': from the samples table (start: wasp),
    else this run's split bam."""
    if starts_at_wasp(config):
        row = wasp_row(wildcards, modality)
        return {'bam': row['bam'], 'bai': row['bai']}
    sb = config['modules']['split_bams']
    donor = wasp_donor(wildcards, modality)
    bam = f"{sb['split_bams_dir']}{split_bams_fold(config, modality)}{donor}.bam".format(
        pool=wildcards.pool, donor=donor)
    return {'bam': bam, 'bai': f"{bam}.bai"}


def wasp_genotype_files(wildcards, modality: str) -> dict:
    return genotype_files(config, wasp_row(wildcards, modality)['genotype_set'])


def wasp_samples_param(modality: str):
    """params of a donor job: the donor (as split), the genotype IDs,
    space- and ';'-separated, the work dir and paired-end."""
    return {
        'donor': lambda wc: wasp_donor(wc, modality),
        'samples': lambda wc: " ".join(wasp_row(wc, modality)['samples']),
        'sample_ids': lambda wc: ";".join(wasp_row(wc, modality)['samples']),
        'work': lambda wc: work_dir(config, modality).format(pool=wc.pool, final_donor=wc.final_donor),
        'paired': "--paired" if paired_end(config, modality) else "",
    }


_WASP_NOTED = set()


def wasp_donors(modality: str) -> list:
    """(pool, donor, final donor) of every WASP job of 'modality': the rows
    of the samples table (start: wasp), else the split donors that are in
    the table."""
    table = samples_table(config)
    if starts_at_wasp(config):
        missing = missing_bams(table)
        if missing:
            raise ValueError(
                "'workflow: start: wasp' filters existing bams; not there:\n  "
                + "\n  ".join(missing))
        return sorted((p, d, r['final_donor']) for (p, d, m), r in table.items() if m == modality)
    pools = split_pools() if EXTERNAL else POOLS
    out, skipped = [], []
    for pool in pools:
        # Every split modality has the donors of the same labels
        donors = get_real_donors({'pool': pool, 'modality': SPLIT_MODALITY})['donor']
        kept, missing = select_donors(table, pool, donors, modality)
        out += [(pool, d, table_row(table, pool, d, modality)['final_donor']) for d in kept]
        skipped += [f"{pool} {d}" for d in missing]
    if skipped and modality not in _WASP_NOTED:
        _WASP_NOTED.add(modality)
        print(f"Warning: split donor(s) not in the WASP samples table "
              f"{wasp_config(config)['samples']}, their {modality} bams not "
              "WASP-filtered: " + ", ".join(skipped), file=sys.stderr)
    if not out:
        raise ValueError(
            f"No split donor is in the WASP samples table {wasp_config(config)['samples']} "
            f"for {modality}: check its pool, donor and modality columns")
    return out


def wasp_donor_stats(wildcards) -> list:
    stats = donor_paths(config, wildcards.modality)['stats']
    return [stats.format(pool=p, final_donor=f) for p, d, f in wasp_donors(wildcards.modality)]


# Genotype sets given as per-chromosome VCFs: WASP's h5 files, once per set.
# snp2h5 reads each VCF as one chromosome (the CHROM of its first record) and
# keeps all its samples and sites; a chromosome of chrom_sizes with no VCF
# is left out (lib/wasp.py: vcf_chrom_files)
rule wasp_snp2h5:
    input:
        vcf = lambda wc: vcf_files(config, wc.gset),
        chrom_sizes = lambda wc: wasp_config(config)['genotypes'][wc.gset]['chrom_sizes'],

    output:
        **{k: f"{converted_dir(config, '{gset}')}{f}" for k, f in H5_FILES.items()}

    # Only the sets given as VCFs are made here
    wildcard_constraints:
        gset = "|".join(map(re.escape, vcf_sets(config))) or r"(?!)",

    log:
        f"{converted_dir(config, '{gset}')}snp2h5.log"

    resources:
        **get_resources("wasp_snp2h5")

    envmodules:
        *WASP_MODULES,

    shell:
        r"""
        set -euo pipefail
        : "${{WASP_HOME:?WASP_HOME is not set: load the WASP module (references: tools: wasp)}}"
        $WASP_HOME/snp2h5/snp2h5 --chrom {input.chrom_sizes} --format vcf \
            --snp_index {output.snp_index} --snp_tab {output.snp_tab} \
            --haplotype {output.haps} {input.vcf} &> {log}
        """


# The samples table's IDs are in their sets' haps.h5, once per run --------
rule wasp_check_inputs:
    input:
        table = wasp_config(config).get('samples') or [],
        names = WASP_NAMES or [],
        h5 = [f for s in genotype_sets(config) for f in genotype_files(config, s).values()],

    output:
        checked_inputs(config)

    params:
        sets = " ".join(
            f"--set {s}={genotype_files(config, s)['haps']},{genotype_files(config, s)['snp_tab']}"
            for s in genotype_sets(config)),
        names = WASP_NAMES_ARG,
        need_bam = "--need_bam" if starts_at_wasp(config) else "",

    resources:
        **get_resources("wasp_check_inputs")

    envmodules:
        *WASP_MODULES,

    shell:
        r"""
        python helper_py_scripts/wasp_check_inputs.py samples --table {input.table} \
            {params.sets} {params.names} {params.need_bam} -o {output}
        """


# The stats of a modality, at the end of the run: one row per donor bam, a
# summary and plots. Its inputs are this run's donors (so it waits for them),
# but it joins every donor's stats file in the modality's folder: donors run
# as batches (one samples table, then another, into the same wasp_dir) all
# end up in it, a later batch doesn't overwrite it with only its own
rule wasp_stats:
    input:
        wasp_donor_stats

    output:
        **run_outputs(config, '{modality}')

    params:
        folder = lambda wc: modality_dir(config, wc.modality),

    resources:
        **get_resources("wasp_stats")

    conda: "../../envs/basic_sctools.yaml"

    shell:
        r"""
        python3 helper_py_scripts/wasp_stats.py summarize {input} --collect {params.folder} \
            --stats {output.stats} --summary {output.summary} --plots {output.plots}
        """
