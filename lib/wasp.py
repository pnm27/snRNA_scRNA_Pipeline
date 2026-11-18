#!/usr/bin/env python3
"""WASP filtering of per-donor bams (stage 'wasp', rules/wasp/).

Each per-donor bam's reads over the donor's heterozygous SNPs are remapped
with the alleles swapped, and kept only if they map back to the same place
(github.com/bmvdgeijn/WASP). What depends on the modality (the remapper,
paired-end reads, duplicates) is in lib/wasp_modality.py:

  cDNA   sc/snRNA-seq and the RNA (GEX) part of multiome: STAR, duplicates
         kept (rules/wasp/rna.smk)
  ATAC   the ATAC part of multiome (or scATAC): bwa-mem, paired-end,
         duplicates marked by cellranger-arc removed (rules/wasp/atac.smk)

This module and rules/wasp/common.smk are shared by both.

The bams are this run's split bams ('workflow: until: wasp'; of the
modalities of 'split_bams: modalities') or the ones listed in the samples
table ('workflow: start: wasp'). The samples table ('modules: wasp:
samples') gives the genotype sample ID(s) of each donor
(find_intersecting_snps.py's --samples). Columns:

  pool          pool (no '/')
  donor         donor (no '/'); names the output bam
  modality      optional: cDNA or ATAC. Empty (or no column): with start:
                wasp, cDNA; after split_bams, every split modality of the
                donor (a row of the donor's modality wins)
  sample        optional: genotype sample ID(s) in the genotype set,
                separated by ';' (or ',' in a tab-separated table); no
                whitespace in an ID. They're written one per line to the
                file given as --samples (WASP takes a file or a
                comma-separated string, and reads the first word of each
                line). Empty or no column: from the sample names table,
                else the donor name
  bam           the per-donor bam (start: wasp only)
  bai           optional; default <bam>.bai
  genotype_set  which 'modules: wasp: genotypes' set (needed with more
                than one set)

The sample names table ('modules: wasp: sample_names', optional) gives
the genotype IDs of pools or donors whose names aren't the ones in the
genotypes, e.g. one donor per pool, named after the pool. Columns:

  pool          pool
  donor         optional: the donor; empty (or no column): every donor
                of the pool
  sample        genotype sample ID(s), as in the samples table

A (pool, donor) row wins over a pool's row.

Tab-separated, or comma-separated for .csv.
"""

import csv
import os
import re

from .modality import normalize_modality, split_modalities


# Separators of the sample IDs in the 'sample' column
SAMPLE_SEP = re.compile(r'[;,]')
_TABLES: dict = {}


def split_samples(text: str) -> list:
    """Sample IDs of a 'sample' cell."""
    return [s.strip() for s in SAMPLE_SEP.split(text or '') if s.strip()]


def _bad_ids(samples: list):
    """Error text of IDs with whitespace, else None."""
    spaced = [s for s in samples if re.search(r'\s', s)]
    if spaced:
        return (f"sample IDs can't contain whitespace ({', '.join(spaced)}): "
                "find_intersecting_snps.py reads the first word of each line of "
                "--samples; separate several IDs with ';' or ','")
    return None


def _read_rows(path: str) -> list:
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter=',' if path.endswith('.csv') else '\t'))


def read_sample_names(path: str) -> dict:
    """{(pool, donor or ''): [sample IDs]} of a sample names table
    (columns pool, sample and optionally donor; '' = every donor of the pool)."""
    rows = _read_rows(path)
    header = list(rows[0].keys()) if rows else []
    missing = [c for c in ('pool', 'sample') if c not in header]
    if not rows or missing:
        raise ValueError(
            f"WASP sample names table {path}: needs a header with the columns "
            "pool, sample (optionally donor), and one row per pool or donor"
            + (f"; missing: {', '.join(missing)}" if rows else ""))
    names, dups, errors = {}, set(), []
    for n, row in enumerate(rows, start=2):
        get = lambda c: (row.get(c) or '').strip()
        pool, donor, samples = get('pool'), get('donor'), split_samples(get('sample'))
        if not pool or not samples:
            errors.append(f"line {n}: empty pool or sample")
            continue
        bad = _bad_ids(samples)
        if bad:
            errors.append(f"line {n}: {bad}")
            continue
        if (pool, donor) in names:
            dups.add(f"{pool} {donor}".strip())
        names[(pool, donor)] = samples
    if dups:
        errors.append(f"listed twice: {', '.join(sorted(dups))}")
    if errors:
        raise ValueError(f"WASP sample names table {path}:\n  " + "\n  ".join(errors))
    return names


def wasp_config(config) -> dict:
    return (config.get('modules') or {}).get('wasp') or {}


def starts_at_wasp(config) -> bool:
    workflow = config.get('workflow')
    return isinstance(workflow, dict) and workflow.get('start') == 'wasp'


def genotype_sets(config) -> list:
    return list((wasp_config(config).get('genotypes') or {}).keys())


def read_samples_table(path: str, need_bam: bool = False, sets=None,
                       names=None) -> dict:
    """{(pool, donor, modality): {'samples': [...], 'bam', 'bai', 'genotype_set'}}.

    modality: 'cDNA', 'ATAC', or '' for every modality (see table_row).
    need_bam: the bam column is required (start: wasp); an empty modality
    is then cDNA.
    sets: the genotype set names; a row without genotype_set takes the only
    one, and a set not in 'sets' is an error.
    names: the sample names table (read_sample_names), for the rows with
    no sample: its (pool, donor) row, else its pool's row, else the donor.
    """
    sets = list(sets or [])
    key = (path, need_bam, tuple(sets), names)
    if key in _TABLES:
        return _TABLES[key]
    name_map = read_sample_names(names) if names else {}
    rows = _read_rows(path)
    header = list(rows[0].keys()) if rows else []
    required = ['pool', 'donor'] + (['bam'] if need_bam else [])
    missing = [c for c in required if c not in header]
    if not rows or missing:
        raise ValueError(
            f"WASP samples table {path}: needs a header with the columns "
            f"{', '.join(required)} (optionally "
            f"{', '.join(c for c in ('modality', 'sample', 'bam', 'bai', 'genotype_set') if c not in required)}),"
            " and one row per donor bam"
            + (f"; missing: {', '.join(missing)}" if rows else ""))
    table, dups, errors = {}, set(), []
    for n, row in enumerate(rows, start=2):
        get = lambda c: (row.get(c) or '').strip()
        pool, donor = get('pool'), get('donor')
        samples = (split_samples(get('sample')) or name_map.get((pool, donor))
                   or name_map.get((pool, '')) or [donor])
        if not pool or not donor:
            errors.append(f"line {n}: empty pool or donor")
            continue
        bad = _bad_ids(samples)
        if bad:
            errors.append(f"line {n}: {bad}")
            continue
        if '/' in pool or '/' in donor:
            errors.append(f"line {n}: pool and donor names can't contain '/' ({pool}, {donor})")
            continue
        modality = get('modality')
        if modality:
            try:
                modality = normalize_modality(modality)
            except ValueError as e:
                errors.append(f"line {n}: {e}")
                continue
        elif need_bam:
            modality = 'cDNA'
        name = f"{pool} {donor}" + (f" {modality}" if modality else "")
        bam = get('bam') or None
        if need_bam and not bam:
            errors.append(f"line {n}: no bam for {name}")
            continue
        gset = get('genotype_set') or None
        if sets:
            if gset is None:
                if len(sets) > 1:
                    errors.append(
                        f"line {n}: no genotype_set for {name}, with "
                        f"several genotype sets ({', '.join(sets)})")
                    continue
                gset = sets[0]
            elif gset not in sets:
                errors.append(
                    f"line {n}: genotype_set '{gset}' isn't one of "
                    f"'modules: wasp: genotypes' ({', '.join(sets)})")
                continue
        if (pool, donor, modality) in table:
            dups.add(name)
        table[(pool, donor, modality)] = {
            'samples': samples,
            'bam': bam,
            'bai': (get('bai') or (f"{bam}.bai" if bam else None)),
            'genotype_set': gset,
        }
    if dups:
        errors.append(f"pool and donor (and modality) listed twice: {', '.join(sorted(dups))}")
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
                              sets=genotype_sets(config), names=sample_names(config))


def sample_names(config):
    """The sample names table ('modules: wasp: sample_names'), or None."""
    return wasp_config(config).get('sample_names') or None


def table_row(table: dict, pool: str, donor: str, modality: str):
    """The row of a donor's bam of 'modality': its own, else the donor's
    row without a modality; None if neither."""
    return table.get((pool, donor, modality)) or table.get((pool, donor, ''))


def row_name(key: tuple) -> str:
    """'<pool> <donor>[ <modality>]' of a samples table key."""
    return " ".join(k for k in key if k)


def table_modalities(table: dict) -> list:
    """Modalities of the rows of a start: wasp samples table."""
    return [m for m in ('cDNA', 'ATAC') if any(k[2] == m for k in table)]


def missing_bams(table: dict) -> list:
    """'<pool> <donor> <modality>: no bam|bai <path>' of the rows whose files aren't there."""
    out = []
    for key, row in sorted(table.items()):
        for k in ('bam', 'bai'):
            if row[k] and not os.path.isfile(row[k]):
                out.append(f"{row_name(key)}: no {k} {row[k]}")
    return out


def select_donors(table: dict, pool: str, donors: list, modality: str = 'cDNA') -> tuple:
    """(donors of 'pool' with a row for 'modality', donors of 'pool' without)."""
    kept = [d for d in donors if table_row(table, pool, d, modality) is not None]
    return kept, [d for d in donors if d not in kept]


def wasp_modalities(config) -> list:
    """Modalities WASP-filtered in this run: those of the samples table
    (start: wasp), else the split modalities ('split_bams: modalities')
    that a row of the samples table applies to."""
    table = samples_table(config)
    if starts_at_wasp(config):
        return table_modalities(table)
    return [m for m in split_modalities(config)
            if any(k[2] in (m, '') for k in table)]


def set_samples(table: dict, gset: str) -> list:
    """The genotype sample IDs of a genotype set, sorted, once each."""
    return sorted({s for row in table.values() if row['genotype_set'] == gset
                   for s in row['samples']})


# Paths --------------------------------------------------------------------

def _dir(config) -> str:
    d = wasp_config(config).get('wasp_dir') or 'wasp_bams/'
    return d if d.endswith('/') else f"{d}/"


def donor_paths(config, modality: str) -> dict:
    """Outputs of a donor's WASP job of 'modality' (templates with {pool}
    and {donor}). ATAC also has the flagstat after WASP, before the
    duplicates are removed."""
    d = f"{_dir(config)}{modality}/{{pool}}/"
    out = {
        'bam': f"{d}{{donor}}.after_wasp.bam",
        'bai': f"{d}{{donor}}.after_wasp.bam.bai",
        'stats': f"{d}stats/{{donor}}.wasp_stats.tsv",
        'flagstat_before': f"{d}stats/{{donor}}.before_wasp.flagstat",
        'flagstat_after': f"{d}stats/{{donor}}.after_wasp.flagstat",
    }
    if modality == 'ATAC':
        out['flagstat_wasp'] = f"{d}stats/{{donor}}.wasp_with_duplicates.flagstat"
    return out


# Remapper of each modality (names its log)
REMAPPER = {'cDNA': 'star', 'ATAC': 'bwa'}


def donor_logs(config, modality: str) -> dict:
    d = f"{_dir(config)}{modality}/{{pool}}/logs/"
    return {
        'intersect': f"{d}{{donor}}.find_intersecting_snps.log",
        'remap': f"{d}{{donor}}.{REMAPPER[modality]}_remap.log",
        'filter': f"{d}{{donor}}.filter_remapped_reads.log",
    }


def work_dir(config, modality: str) -> str:
    """Per-donor scratch dir (template), removed when the job ends."""
    d = wasp_config(config).get('work_dir') or 'wasp_tmp/'
    return f"{d if d.endswith('/') else d + '/'}{modality}/{{pool}}/{{donor}}/"


def run_outputs(config, modality: str) -> dict:
    """The run's stats of 'modality' (or '{modality}', a wildcard)."""
    d = f"{_dir(config)}{modality}/"
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
