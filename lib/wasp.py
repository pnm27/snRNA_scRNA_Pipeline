#!/usr/bin/env python3
"""WASP filtering of per-donor bams (stage 'wasp', rules/wasp.smk).

For sc/snRNA-seq and the RNA (GEX) part of multiome. Each per-donor bam's
reads over the donor's heterozygous SNPs are remapped with the alleles
swapped (STAR), and kept only if they map back to the same place
(github.com/bmvdgeijn/WASP). Duplicates are kept.

The bams are this run's split bams ('workflow: until: wasp') or the ones
listed in the samples table ('workflow: start: wasp'). The samples table
('modules: wasp: samples') gives the genotype sample ID(s) of each donor
(find_intersecting_snps.py's --samples). Columns:

  pool          pool (no '/')
  donor         donor (no '/'); names the output bam
  sample        genotype sample ID(s) in the genotype set, separated by
                ';' (or ',' in a tab-separated table); no whitespace in
                an ID. They're written one per line to the file given
                as --samples (WASP takes a file or a comma-separated
                string, and reads the first word of each line)
  bam           the per-donor bam (start: wasp only)
  bai           optional; default <bam>.bai
  genotype_set  which 'modules: wasp: genotypes' set (needed with more
                than one set)

Tab-separated, or comma-separated for .csv.
"""

import csv
import os
import re


# Separators of the sample IDs in the 'sample' column
SAMPLE_SEP = re.compile(r'[;,]')
_TABLES: dict = {}


def wasp_config(config) -> dict:
    return (config.get('modules') or {}).get('wasp') or {}


def starts_at_wasp(config) -> bool:
    workflow = config.get('workflow')
    return isinstance(workflow, dict) and workflow.get('start') == 'wasp'


def genotype_sets(config) -> list:
    return list((wasp_config(config).get('genotypes') or {}).keys())


def read_samples_table(path: str, need_bam: bool = False, sets=None) -> dict:
    """{(pool, donor): {'samples': [...], 'bam', 'bai', 'genotype_set'}}.

    need_bam: the bam column is required (start: wasp).
    sets: the genotype set names; a row without genotype_set takes the only
    one, and a set not in 'sets' is an error.
    """
    sets = list(sets or [])
    key = (path, need_bam, tuple(sets))
    if key in _TABLES:
        return _TABLES[key]
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f, delimiter=',' if path.endswith('.csv') else '\t'))
    header = list(rows[0].keys()) if rows else []
    required = ['pool', 'donor', 'sample'] + (['bam'] if need_bam else [])
    missing = [c for c in required if c not in header]
    if not rows or missing:
        raise ValueError(
            f"WASP samples table {path}: needs a header with the columns "
            f"{', '.join(required)} (optionally "
            f"{', '.join(c for c in ('bam', 'bai', 'genotype_set') if c not in required)}),"
            " and one row per donor bam"
            + (f"; missing: {', '.join(missing)}" if rows else ""))
    table, dups, errors = {}, set(), []
    for n, row in enumerate(rows, start=2):
        get = lambda c: (row.get(c) or '').strip()
        pool, donor = get('pool'), get('donor')
        samples = [s.strip() for s in SAMPLE_SEP.split(get('sample')) if s.strip()]
        if not pool or not donor or not samples:
            errors.append(f"line {n}: empty pool, donor or sample")
            continue
        spaced = [s for s in samples if re.search(r'\s', s)]
        if spaced:
            errors.append(
                f"line {n}: sample IDs can't contain whitespace ({', '.join(spaced)}): "
                "find_intersecting_snps.py reads the first word of each line of "
                "--samples; separate several IDs with ';' or ','")
            continue
        if '/' in pool or '/' in donor:
            errors.append(f"line {n}: pool and donor names can't contain '/' ({pool}, {donor})")
            continue
        bam = get('bam') or None
        if need_bam and not bam:
            errors.append(f"line {n}: no bam for {pool} {donor}")
            continue
        gset = get('genotype_set') or None
        if sets:
            if gset is None:
                if len(sets) > 1:
                    errors.append(
                        f"line {n}: no genotype_set for {pool} {donor}, with "
                        f"several genotype sets ({', '.join(sets)})")
                    continue
                gset = sets[0]
            elif gset not in sets:
                errors.append(
                    f"line {n}: genotype_set '{gset}' isn't one of "
                    f"'modules: wasp: genotypes' ({', '.join(sets)})")
                continue
        if (pool, donor) in table:
            dups.add(f"{pool} {donor}")
        table[(pool, donor)] = {
            'samples': samples,
            'bam': bam,
            'bai': (get('bai') or (f"{bam}.bai" if bam else None)),
            'genotype_set': gset,
        }
    if dups:
        errors.append(f"pool and donor listed twice: {', '.join(sorted(dups))}")
    if errors:
        raise ValueError(f"WASP samples table {path}:\n  " + "\n  ".join(errors))
    _TABLES[key] = table
    return table


def samples_table(config) -> dict:
    """The run's samples table (see read_samples_table)."""
    path = wasp_config(config).get('samples')
    if not path:
        raise ValueError("The wasp stage needs a samples table: set 'modules: wasp: samples'")
    return read_samples_table(path, need_bam=starts_at_wasp(config),
                              sets=genotype_sets(config))


def missing_bams(table: dict) -> list:
    """'<pool> <donor>: no bam|bai <path>' of the rows whose files aren't there."""
    out = []
    for (pool, donor), row in sorted(table.items()):
        for k in ('bam', 'bai'):
            if row[k] and not os.path.isfile(row[k]):
                out.append(f"{pool} {donor}: no {k} {row[k]}")
    return out


def select_donors(table: dict, pool: str, donors: list) -> tuple:
    """(donors of 'pool' in the table, donors of 'pool' not in it)."""
    kept = [d for d in donors if (pool, d) in table]
    return kept, [d for d in donors if (pool, d) not in table]


def set_samples(table: dict, gset: str) -> list:
    """The genotype sample IDs of a genotype set, sorted, once each."""
    return sorted({s for row in table.values() if row['genotype_set'] == gset
                   for s in row['samples']})


# Paths --------------------------------------------------------------------

def _dir(config) -> str:
    d = wasp_config(config).get('wasp_dir') or 'wasp_bams/'
    return d if d.endswith('/') else f"{d}/"


def donor_paths(config) -> dict:
    """Outputs of a donor's WASP job (templates with {pool} and {donor})."""
    d = f"{_dir(config)}{{pool}}/"
    return {
        'bam': f"{d}{{donor}}.after_wasp.bam",
        'bai': f"{d}{{donor}}.after_wasp.bam.bai",
        'stats': f"{d}stats/{{donor}}.wasp_stats.tsv",
        'flagstat_before': f"{d}stats/{{donor}}.before_wasp.flagstat",
        'flagstat_after': f"{d}stats/{{donor}}.after_wasp.flagstat",
    }


def donor_logs(config) -> dict:
    d = f"{_dir(config)}{{pool}}/logs/"
    return {
        'intersect': f"{d}{{donor}}.find_intersecting_snps.log",
        'star': f"{d}{{donor}}.star_remap.log",
        'filter': f"{d}{{donor}}.filter_remapped_reads.log",
    }


def work_dir(config) -> str:
    """Per-donor scratch dir (template), removed when the job ends."""
    d = wasp_config(config).get('work_dir') or 'wasp_tmp/'
    return f"{d if d.endswith('/') else d + '/'}{{pool}}/{{donor}}/"


def run_outputs(config) -> dict:
    d = _dir(config)
    return {
        'stats': f"{d}wasp_stats.tsv",
        'summary': f"{d}wasp_summary.tsv",
        'plots': f"{d}wasp_plots.pdf",
    }


def checked_inputs(config) -> str:
    return f"{_dir(config)}checked_inputs.txt"


H5_FILES = {'snp_index': 'snp_index.h5', 'snp_tab': 'snp_tab.h5', 'haps': 'haps.h5'}


def converted_dir(config, gset: str) -> str:
    """Where a set given as VCF/BCF is converted to WASP's h5 files."""
    return f"{_dir(config)}genotypes/{gset}/"


def genotype_files(config, gset: str) -> dict:
    """snp_index, snp_tab and haps h5 files of a genotype set."""
    g = (wasp_config(config).get('genotypes') or {}).get(gset)
    if g is None:
        raise ValueError(f"No genotype set '{gset}' in 'modules: wasp: genotypes'")
    d = g.get('h5_dir') or converted_dir(config, gset)
    d = d if d.endswith('/') else f"{d}/"
    return {k: f"{d}{f}" for k, f in H5_FILES.items()}


def vcf_sets(config) -> list:
    """Genotype sets given as VCF/BCF (converted with snp2h5)."""
    return [s for s, g in (wasp_config(config).get('genotypes') or {}).items()
            if not (g or {}).get('h5_dir')]


def vcf_files(config, gset: str) -> list:
    vcf = wasp_config(config)['genotypes'][gset].get('vcf')
    return [vcf] if isinstance(vcf, str) else list(vcf or [])


def remap_reference(config) -> tuple:
    """(STAR genome dir, gtf or None) to remap with.

    'modules: wasp: remap: genome' (and gtf), else the run's aligner
    reference: STARsolo's genome_pick (genome and gtf), or the 'star' index
    of a cellranger / cellranger-arc reference (its annotations are in the
    index). A run that starts at wasp names the genome.
    """
    r = wasp_config(config).get('remap') or {}
    if r.get('genome'):
        return r['genome'], r.get('gtf')
    if starts_at_wasp(config):
        raise ValueError(
            "'workflow: start: wasp' needs 'modules: wasp: remap: genome' "
            "(the STAR index of the genome the bams were aligned to)")
    aligner = str((config.get('workflow') or {}).get('aligner')).lower()
    mods = config['modules']
    if aligner == 'starsolo':
        pick = mods['STARsolo']['genome_pick']
        return pick['genome'], r.get('gtf') or pick.get('gtf')
    section = {'cellranger': 'cellranger', 'cellranger_arc': 'cellranger_arc'}.get(aligner)
    if section is None:
        raise ValueError(f"Unrecognised aligner: {aligner}")
    genome = mods[section]['ref']['genome']
    return f"{genome.rstrip('/')}/star", r.get('gtf')
