# WASP filtering of per-donor bams ('workflow: until: wasp', or 'start:
# wasp' for bams listed in the samples table): what every modality shares.
# See lib/wasp.py and docs/source/sub_workflows.md.
#
#   common.smk  genotypes (snp2h5), input checks, the stats (a run of
#               their own: target wasp_stats_all), the
#               per-donor prepare and finish jobs, and their input functions
#   rna.smk     cDNA (sc/snRNA-seq, multiome GEX): STAR remapping
#   atac.smk    ATAC (multiome ATAC, scATAC): bwa-mem remapping
#
# A pipeline of one modality includes common.smk and its own file. Per donor
# bam, as github.com/bmvdgeijn/WASP's mapping pipeline:
# find_intersecting_snps.py (reads over the donor's SNPs, the --samples IDs
# of the samples table), remapping of those reads with the alleles swapped,
# filter_remapped_reads.py (keeps the reads that map back to the same
# place), merged with the reads over no SNP. Three jobs per donor bam, so a
# failed one is rerun alone (the files the next one reads are kept):
#
#   wasp_prepare_<rna|atac>  wasp_donor.sh prepare: chromosome names,
#                            find_intersecting_snps.py
#   wasp_remap_<rna|atac>    the remapping (rna.smk, atac.smk)
#   wasp_finish_<rna|atac>   wasp_donor.sh finish: filter_remapped_reads.py,
#                            merge, sort, duplicates (ATAC), index, flagstat,
#                            stats; removes the donor's work dir
#
# The files between them are temp() in the donor's work dir ('work_dir').
# WASP comes from its env module ('references: tools: wasp'), whose scripts
# are in $WASP_HOME.
#
# At most 'wasp_slots' jobs at a time (each holds one): the temporary bams
# are large, e.g. --resources wasp_slots=100.

from lib.wasp import (
    samples_table, starts_at_wasp, select_donors, missing_bams, table_row, donor_of,
    donor_paths, donor_logs, work_dir, modality_dir, checked_inputs,
    intersect_paths, remap_fastqs, remapped_bam, remap_scratch,
    genotype_files, genotype_sets, vcf_sets, vcf_files, converted_dir,
    wasp_config, sample_names, H5_FILES,
)
from lib.wasp_modality import paired_end, drop_duplicates
import glob as _glob
import os
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


# A donor's stats file in a modality's folder (as wasp_stats.py --collect)
DONOR_STATS_GLOB = os.path.join('*', 'stats', '*.wasp_stats.tsv')


def wasp_stats_modalities() -> list:
    """The modalities of wasp_dir that have donors' stats files."""
    out = [m for m in WASP_RULE_SUFFIX
           if _glob.glob(os.path.join(modality_dir(config, m), DONOR_STATS_GLOB))]
    if not out:
        raise ValueError(f"No donor's WASP stats in any modality's folder of "
                         f"{modality_dir(config, '')} ({DONOR_STATS_GLOB}): run the donors first")
    return out


# The per-donor jobs before and after the remapping, one rule of each per
# modality ---------------------------------------------------------------
WASP_RULE_SUFFIX = {'cDNA': 'rna', 'ATAC': 'atac'}

for _modality, _suffix in WASP_RULE_SUFFIX.items():
    _paired = paired_end(config, _modality)

    # Reads over the donor's SNPs (to remap) and the others
    rule:
        name: f"wasp_prepare_{_suffix}"

        input:
            unpack(partial(wasp_input_bam, modality=_modality)),
            unpack(partial(wasp_genotype_files, modality=_modality)),
            checked = checked_inputs(config),

        output:
            **{k: temp(p) for k, p in intersect_paths(config, _modality, _paired).items()}

        params:
            **wasp_samples_param(_modality),

        log:
            donor_logs(config, _modality)['intersect']

        resources:
            **get_resources(f"wasp_prepare_{_suffix}"),
            wasp_slots = 1,

        envmodules:
            *WASP_MODULES,
            "samtools",

        shell:
            r"""
            bash helper_sh_scripts/wasp_donor.sh prepare --work {params.work} --bam {input.bam} \
                --samples {params.samples:q} --snp_index {input.snp_index} \
                --snp_tab {input.snp_tab} --haps {input.haps} {params.paired} \
                --log {log}
            """

    # The remapped reads that map back to the same place, with the others
    rule:
        name: f"wasp_finish_{_suffix}"

        input:
            unpack(partial(wasp_input_bam, modality=_modality)),
            to_remap = intersect_paths(config, _modality, _paired)['to_remap'],
            keep = intersect_paths(config, _modality, _paired)['keep'],
            remapped = remapped_bam(config, _modality),

        output:
            **donor_paths(config, _modality)

        params:
            **wasp_samples_param(_modality),
            # the prepare job's log, which the stats read
            intersect_log = (lambda m: lambda wc: donor_logs(config, m)['intersect'].format(
                pool=wc.pool, final_donor=wc.final_donor))(_modality),
            # ATAC: the flagstat before the duplicates are removed, and removing them
            flagstat_wasp = lambda wc, output: (
                f"--flagstat_wasp {output.flagstat_wasp}" if 'flagstat_wasp' in output.keys() else ""),
            drop_duplicates = (
                "--drop_duplicates" if _modality == 'ATAC' and drop_duplicates(config) else ""),
            modality = _modality,
            threads = lambda wc, resources: resources.cpus_per_task,

        log:
            donor_logs(config, _modality)['filter']

        resources:
            **get_resources(f"wasp_finish_{_suffix}"),
            wasp_slots = 1,

        envmodules:
            *WASP_MODULES,
            "samtools",

        shell:
            r"""
            bash helper_sh_scripts/wasp_donor.sh finish --work {params.work} --bam {input.bam} \
                --remapped {input.remapped} --out {output.bam} \
                --threads {params.threads} --filter_log {log} \
                --intersect_log {params.intersect_log} --flagstat_before {output.flagstat_before} \
                --flagstat_after {output.flagstat_after} {params.flagstat_wasp} \
                {params.drop_duplicates} {params.paired} \
                --pool {wildcards.pool} --donor {params.donor} --final_donor {wildcards.final_donor} \
                --sample {params.sample_ids:q} --modality {params.modality} --stats {output.stats}
            """


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


# The stats of each modality: one row per donor bam, a summary and plots
# (<wasp_dir><modality>/wasp_stats.tsv, wasp_summary.tsv, wasp_plots.pdf), of
# every donor's stats file in the modality's folder (all batches run into
# the wasp_dir). Not a target of the donors' runs, so a batch neither remakes
# them nor waits for them: a run of their own, after (or between) the
# batches, with any batch's config (only its wasp_dir is used):
#
#   snakemake --configfile <config> wasp_stats --use-conda ...
#
# No files declared: the donors' stats files would tie it to their jobs
# (whose donors another batch's samples table doesn't have), and it runs
# every time it's asked for. Each file is written to a temporary name and
# moved into place. Each run has its own log, named by its start time
# (<wasp_dir>logs/wasp_stats.<YYYYmmdd_HHMMSS>.log): when, the pipeline's
# commit and config, and each modality's stats files used (with their
# times) and files written
rule wasp_stats:
    params:
        folders = lambda wc: " ".join(modality_dir(config, m) for m in wasp_stats_modalities()),
        log_dir = f"{modality_dir(config, '')}logs/",
        configfiles = " ".join(os.path.abspath(f) for f in workflow.overwrite_configfiles) or "none",

    resources:
        **get_resources("wasp_stats")

    conda: "../../envs/basic_sctools.yaml"

    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.log_dir}
        log={params.log_dir}wasp_stats.$(date +%Y%m%d_%H%M%S).log
        echo "wasp_stats: log $log" >&2
        exec 3>&2 >>"$log" 2>&1
        trap 'echo "wasp_stats: failed (exit $?), see $log" >&3' ERR
        echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] wasp_stats started on $(hostname)"
        echo "pipeline: $PWD, commit $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
        echo "config: {params.configfiles}"
        for folder in {params.folders}; do
            echo "== $folder"
            python3 helper_py_scripts/wasp_stats.py summarize --collect "$folder" -o "$folder"                 --log "$log"
        done
        echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] wasp_stats done"
        """
