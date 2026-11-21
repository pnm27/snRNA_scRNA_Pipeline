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
  donor         donor (no '/') as the final demultiplexing named it (the
                split bams' name)
  final_donor   optional: the donor after swap correction (no '/'); names
                the output bam. Empty (or no column): from the sample
                names table, else the donor (no swap corrections; a
                warning says so). Once per pool and modality
  modality      optional: cDNA or ATAC. Empty (or no column): with start:
                wasp, cDNA; after split_bams, every split modality of the
                donor (a row of the donor's modality wins)
  vcfID         optional: genotype sample ID(s) in the genotype set (of
                the final donor; 'sample' is read too),
                separated by ';' (or ',' in a tab-separated table); no
                whitespace in an ID. They're written one per line to the
                file given as --samples (WASP takes a file or a
                comma-separated string, and reads the first word of each
                line). Empty or no column: from the sample names table,
                else the final donor's name
  bam           the per-donor bam (start: wasp only)
  bai           optional; default <bam>.bai
  genotype_set  which 'modules: wasp: genotypes' set (needed with more
                than one set)

The sample names table ('modules: wasp: sample_names', optional) gives
the final donor and genotype IDs of pools or donors whose names aren't
the ones in the genotypes, e.g. one donor per pool, named after the pool,
or donors whose labels were swapped. Columns:

  pool          pool
  donor         optional: the donor, as the final demultiplexing named
                it; empty (or no column): every donor of the pool
  final_donor   optional: the donor after swap correction (the swap
                corrections are only applied to the h5ad files so far);
                empty: the donor. No final_donor column here nor in the
                samples table: a warning, and final_donor = donor (no swap
                corrections assumed)
  vcfID         genotype sample ID(s) of the final donor, as in the
                samples table ('sample' is read too)

A (pool, donor) row wins over a pool's row.

Tab-separated, or comma-separated for .csv.
"""

import csv
import os
import re
import sys

from .modality import normalize_modality, split_modalities


# Separators of the sample IDs in the 'sample' column
SAMPLE_SEP = re.compile(r'[;,]')
# Names of the genotype ID column ('sample': tables from before vcfID)
SAMPLE_COLUMNS = ('vcfID', 'sample')
_TABLES: dict = {}
_WARNED: set = set()


def _warn(message: str):
    """Prints a warning once per process."""
    if message not in _WARNED:
        _WARNED.add(message)
        print(f"Warning: {message}", file=sys.stderr)


def _sample_column(header: list, what: str):
    """The genotype ID column of a table's header, or None."""
    found = [c for c in SAMPLE_COLUMNS if c in header]
    if len(found) > 1:
        raise ValueError(f"{what}: both vcfID and sample columns; keep one (vcfID)")
    return found[0] if found else None


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


def read_sample_names(path: str) -> tuple:
    """({(pool, donor or ''): {'samples': [...], 'final_donor': str}},
    whether the table has a final_donor column) of a sample names table
    (columns pool, vcfID (or sample), optionally donor and final_donor;
    donor '' = every donor of the pool, final_donor '' = the donor)."""
    rows = _read_rows(path)
    header = list(rows[0].keys()) if rows else []
    what = f"WASP sample names table {path}"
    sample_col = _sample_column(header, what)
    missing = [c for c, ok in (('pool', 'pool' in header), ('vcfID', sample_col)) if not ok]
    if not rows or missing:
        raise ValueError(
            f"{what}: needs a header with the columns pool, vcfID (optionally "
            "donor, final_donor), and one row per pool or donor"
            + (f"; missing: {', '.join(missing)}" if rows else ""))
    has_final = 'final_donor' in header
    names, dups, errors = {}, set(), []
    for n, row in enumerate(rows, start=2):
        get = lambda c: (row.get(c) or '').strip()
        pool, donor, samples = get('pool'), get('donor'), split_samples(get(sample_col))
        final = get('final_donor')
        if not pool or not samples:
            errors.append(f"line {n}: empty pool or {sample_col}")
            continue
        bad = _bad_ids(samples)
        if bad:
            errors.append(f"line {n}: {bad}")
            continue
        if '/' in final:
            errors.append(f"line {n}: final_donor can't contain '/' ({final})")
            continue
        if (pool, donor) in names:
            dups.add(f"{pool} {donor}".strip())
        names[(pool, donor)] = {'samples': samples, 'final_donor': final}
    if dups:
        errors.append(f"listed twice: {', '.join(sorted(dups))}")
    if errors:
        raise ValueError(f"{what}:\n  " + "\n  ".join(errors))
    return names, has_final


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
    no vcfID or final_donor: its (pool, donor) row, else its pool's row,
    else the donor.
    Each row also has 'final_donor' (the swap-corrected donor, which names
    the outputs); a pool's final donors of a modality are unique.
    """
    sets = list(sets or [])
    key = (path, need_bam, tuple(sets), names)
    if key in _TABLES:
        return _TABLES[key]
    name_map, names_final = read_sample_names(names) if names else ({}, False)
    rows = _read_rows(path)
    header = list(rows[0].keys()) if rows else []
    what = f"WASP samples table {path}"
    sample_col = _sample_column(header, what) or 'vcfID'
    required = ['pool', 'donor'] + (['bam'] if need_bam else [])
    missing = [c for c in required if c not in header]
    if not rows or missing:
        raise ValueError(
            f"WASP samples table {path}: needs a header with the columns "
            f"{', '.join(required)} (optionally "
            f"{', '.join(c for c in ('final_donor', 'modality', 'vcfID', 'bam', 'bai', 'genotype_set') if c not in required)}),"
            " and one row per donor bam"
            + (f"; missing: {', '.join(missing)}" if rows else ""))
    has_final = 'final_donor' in header
    if not has_final and not names_final:
        _warn(f"{what}" + (f" and sample names table {names} have" if names else " has")
              + " no final_donor column: no swap corrections assumed (final_donor = donor)")
    table, dups, errors, unset = {}, set(), [], []
    for n, row in enumerate(rows, start=2):
        get = lambda c: (row.get(c) or '').strip()
        pool, donor = get('pool'), get('donor')
        named = name_map.get((pool, donor)) or name_map.get((pool, '')) or {}
        samples = split_samples(get(sample_col)) or named.get('samples')
        final = (get('final_donor') or (name_map.get((pool, donor)) or {}).get('final_donor')
                 or (name_map.get((pool, '')) or {}).get('final_donor') or '')
        if not pool or not donor:
            errors.append(f"line {n}: empty pool or donor")
            continue
        if not final:
            final = donor
            if has_final or names_final:
                unset.append(f"{pool} {donor}")
        if '/' in final:
            errors.append(f"line {n}: final_donor can't contain '/' ({final})")
            continue
        samples = samples or [final]
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
            'final_donor': final,
            'samples': samples,
            'bam': bam,
            'bai': (get('bai') or (f"{bam}.bai" if bam else None)),
            'genotype_set': gset,
        }
    if dups:
        errors.append(f"pool and donor (and modality) listed twice: {', '.join(sorted(dups))}")
    if not errors:
        errors += [f"{pool}: final_donor {final} ({modality}) given to several "
                   f"donors ({', '.join(donors)})"
                   for (pool, final, modality), donors in sorted(_final_donors(table).items())
                   if len(donors) > 1]
    if errors:
        raise ValueError(f"{what}:\n  " + "\n  ".join(errors))
    if unset:
        _warn(f"{what}: no final_donor for {', '.join(unset)}: final_donor = donor")
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


MODALITIES = ('cDNA', 'ATAC')


_FINALS: dict = {}


def _final_donors(table: dict) -> dict:
    """{(pool, final_donor, modality): [donors]} of each modality's rows."""
    cached = _FINALS.get(id(table))
    if cached is not None and cached[0] is table:
        return cached[1]
    out = {}
    for pool, donor in sorted({k[:2] for k in table}):
        for m in MODALITIES:
            row = table_row(table, pool, donor, m)
            if row is not None:
                out.setdefault((pool, row['final_donor'], m), []).append(donor)
    _FINALS[id(table)] = (table, out)
    return out


def donor_of(table: dict, pool: str, final_donor: str, modality: str):
    """The donor (split bam name) whose final donor is 'final_donor', or None."""
    donors = _final_donors(table).get((pool, final_donor, modality))
    return donors[0] if donors else None


def row_name(key: tuple) -> str:
    """'<pool> <donor>[ <modality>]' of a samples table key."""
    return " ".join(k for k in key if k)


def table_modalities(table: dict) -> list:
    """Modalities of the rows of a start: wasp samples table."""
    return [m for m in MODALITIES if any(k[2] == m for k in table)]


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
    and {final_donor}, the swap-corrected donor). ATAC also has the flagstat after WASP, before the
    duplicates are removed."""
    d = f"{_dir(config)}{modality}/{{pool}}/"
    out = {
        'bam': f"{d}{{final_donor}}.after_wasp.bam",
        'bai': f"{d}{{final_donor}}.after_wasp.bam.bai",
        'stats': f"{d}stats/{{final_donor}}.wasp_stats.tsv",
        'flagstat_before': f"{d}stats/{{final_donor}}.before_wasp.flagstat",
        'flagstat_after': f"{d}stats/{{final_donor}}.after_wasp.flagstat",
    }
    if modality == 'ATAC':
        out['flagstat_wasp'] = f"{d}stats/{{final_donor}}.wasp_with_duplicates.flagstat"
    return out


# Remapper of each modality (names its log)
REMAPPER = {'cDNA': 'star', 'ATAC': 'bwa'}


def donor_logs(config, modality: str) -> dict:
    d = f"{_dir(config)}{modality}/{{pool}}/logs/"
    return {
        'intersect': f"{d}{{final_donor}}.find_intersecting_snps.log",
        'remap': f"{d}{{final_donor}}.{REMAPPER[modality]}_remap.log",
        'filter': f"{d}{{final_donor}}.filter_remapped_reads.log",
    }


def work_dir(config, modality: str) -> str:
    """Per-donor scratch dir (template): made by the donor's prepare job,
    kept between its jobs, removed when its finish job ends well."""
    d = wasp_config(config).get('work_dir') or 'wasp_tmp/'
    return f"{d if d.endswith('/') else d + '/'}{modality}/{{pool}}/{{final_donor}}/"


def intersect_paths(config, modality: str, paired: bool) -> dict:
    """Outputs of a donor's 'wasp_prepare_*' job (helper_sh_scripts/wasp_donor.sh
    prepare; templates, in the work dir): the reads over SNPs, as a bam (for
    filter_remapped_reads.py) and as the fastqs to remap, and the others."""
    d = f"{work_dir(config, modality)}intersect/"
    out = {'to_remap': f"{d}to.remap.bam", 'keep': f"{d}keep.bam"}
    if paired:
        out.update(fq1=f"{d}remap.fq1.gz", fq2=f"{d}remap.fq2.gz")
    else:
        out['fq'] = f"{d}remap.fq.gz"
    return out


def remap_fastqs(config, modality: str, paired: bool) -> list:
    """The fastqs the remap job reads (fq1 fq2, or fq)."""
    return [p for k, p in intersect_paths(config, modality, paired).items() if k.startswith('fq')]


def remapped_bam(config, modality: str) -> str:
    """Output of a donor's remap job (template, in the work dir)."""
    return f"{work_dir(config, modality)}remapped.bam"


def remap_scratch(config, modality: str) -> str:
    """The aligner's own files of a remap job (template): cleared when it
    starts, removed when it ends well, kept when it fails."""
    return f"{work_dir(config, modality)}remap/"


def modality_dir(config, modality: str) -> str:
    """A modality's folder in wasp_dir: its pools' folders and its stats."""
    return f"{_dir(config)}{modality}/"


def run_outputs(config, modality: str) -> dict:
    """The stats of 'modality' (or '{modality}', a wildcard): of every donor
    in its folder, all batches (rule wasp_stats)."""
    d = modality_dir(config, modality)
    return {
        'stats': f"{d}wasp_stats.tsv",
        'summary': f"{d}wasp_summary.tsv",
        'plots': f"{d}wasp_plots.pdf",
    }


def checked_inputs(config) -> str:
    return f"{_dir(config)}checked_inputs.txt"


H5_FILES = {'snp_index': 'snp_index.h5', 'snp_tab': 'snp_tab.h5', 'haps': 'haps.h5'}


def converted_dir(config, gset: str) -> str:
    """Where a set given as per-chromosome VCFs is converted to WASP's h5 files."""
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
    """Genotype sets given as per-chromosome VCFs (converted with snp2h5)."""
    return [s for s, g in (wasp_config(config).get('genotypes') or {}).items()
            if not (g or {}).get('h5_dir')]


VCF_PATTERN = re.compile(r'.*\{chrom\}.*\.vcf(\.gz)?$')


def read_chrom_sizes(path: str) -> list:
    """Chromosome names (first column) of a chrom sizes file, in order."""
    if not os.path.isfile(path):
        raise ValueError(f"No chrom_sizes file {path}")
    with open(path) as f:
        return [line.split()[0] for line in f if line.strip() and not line.startswith('#')]


def vcf_chrom_files(config, gset: str) -> dict:
    """{chromosome: VCF} of a set's 'vcf' pattern, filled with each chromosome
    of its 'chrom_sizes' whose file exists (in chrom_sizes' order).

    snp2h5 reads one chromosome per file (the CHROM of its first record) and
    only text VCF (.vcf.gz or .vcf), so the set is a pattern with {chrom},
    e.g. myVcf_{chrom}.vcf.gz with chrom_sizes naming chr1, chr2, ...
    Chromosomes without a file (e.g. chrM, alt contigs) are left out, with a
    warning; none at all is an error."""
    g = wasp_config(config)['genotypes'][gset]
    pattern, sizes = g.get('vcf'), g.get('chrom_sizes')
    if not isinstance(pattern, str) or not VCF_PATTERN.fullmatch(pattern):
        raise ValueError(
            f"Genotype set '{gset}': 'vcf' must be a pattern with {{chrom}} of "
            f"per-chromosome VCFs (.vcf.gz or .vcf), e.g. myVcf_{{chrom}}.vcf.gz, "
            f"not {pattern!r}: snp2h5 reads one chromosome per text VCF")
    found, missing = {}, []
    for chrom in read_chrom_sizes(sizes):
        path = pattern.replace('{chrom}', chrom)
        if os.path.exists(path):
            found[chrom] = path
        else:
            missing.append(chrom)
    if not found:
        raise ValueError(
            f"Genotype set '{gset}': no VCF {pattern} for the chromosomes of "
            f"{sizes} (e.g. {pattern.replace('{chrom}', missing[0]) if missing else pattern}); "
            "{chrom} is filled with chrom_sizes' names (chr1 or 1)")
    if missing:
        _warn(f"Genotype set '{gset}': no VCF {pattern} for {', '.join(missing)} "
              f"(in {sizes}): left out of its h5 files")
    return found


def vcf_files(config, gset: str) -> list:
    """The set's per-chromosome VCFs (vcf_chrom_files)."""
    return list(vcf_chrom_files(config, gset).values())
