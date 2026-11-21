#!/usr/bin/env python3
"""What WASP does differently per modality (see lib/wasp.py).

  cDNA  'modules: wasp: rna': STAR remapping (genome, gtf; the STAR module
        'references: tools: wasp_star'), paired_end (false for 10x 3', where
        only R2 has the cDNA); duplicates kept (rules/wasp/rna.smk)
  ATAC  'modules: wasp: atac': bwa-mem remapping (bwa_index, as
        cellranger-arc), always paired-end; with 'drop_duplicates: marked'
        the reads cellranger-arc marked as duplicates (flag 0x400, per cell
        barcode) are removed after WASP (rules/wasp/atac.smk)

A pipeline of one modality (e.g. sc/snRNA-seq only) needs only its own
part of this module, rules/wasp/common.smk and its rules/wasp/<modality>.smk.
"""

from .wasp import wasp_config, starts_at_wasp, _warn


def rna_config(config) -> dict:
    return wasp_config(config).get('rna') or {}


def atac_config(config) -> dict:
    return wasp_config(config).get('atac') or {}


def paired_end(config, modality: str) -> bool:
    """Whether find_intersecting_snps.py reads pairs (--is_paired_end)."""
    if modality == 'ATAC':
        return True
    return bool(rna_config(config).get('paired_end'))


# cDNA (RNA) ---------------------------------------------------------------

def remap_reference(config) -> tuple:
    """(STAR genome dir, gtf or None) to remap the cDNA with.

    'modules: wasp: rna: remap: genome' (and gtf), else the run's aligner
    reference: STARsolo's genome_pick (genome and gtf), or the 'star' index
    of a cellranger / cellranger-arc reference (its annotations are in the
    index). A run that starts at wasp names the genome.
    """
    r = rna_config(config).get('remap') or {}
    if r.get('genome'):
        return r['genome'], r.get('gtf')
    if starts_at_wasp(config):
        raise ValueError(
            "'workflow: start: wasp' needs 'modules: wasp: rna: remap: genome' "
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


def remap_star_module(config) -> str:
    """Env module of the STAR that remaps the cDNA: 'references: tools:
    wasp_star', else the aligner's 'star'.

    A STAR index loads only in a STAR compatible with the one that built it,
    so the remap's STAR is set apart from the aligner's: e.g. the 'star'
    folder of a cellranger / cellranger-arc reference (the default remap
    genome with those aligners) was built by that cellranger's STAR. Using
    'star' for it gives a warning.
    """
    tools = (config.get('references') or {}).get('tools') or {}
    if tools.get('wasp_star'):
        return tools['wasp_star']
    aligner = str((config.get('workflow') or {}).get('aligner')).lower()
    if (aligner in ('cellranger', 'cellranger_arc') and not starts_at_wasp(config)
            and not (rna_config(config).get('remap') or {}).get('genome')):
        _warn(f"'references: tools: wasp_star' isn't set: the WASP cDNA remap loads "
              f"the {aligner} reference's STAR index with 'star' ({tools.get('star')}). "
              "If STAR refuses it (\"Genome version ... is INCOMPATIBLE\"), set wasp_star "
              "to a STAR module that can read that index")
    return tools.get('star')


# ATAC ---------------------------------------------------------------------

DUPLICATE_MODES = ('marked', 'none')


def bwa_index(config) -> str:
    """bwa index (its fasta) to remap the ATAC reads with.

    'modules: wasp: atac: remap: bwa_index', else the cellranger-arc
    reference's fasta/genome.fa (indexed for bwa by cellranger-arc
    mkref). A run that starts at wasp names it.
    """
    r = atac_config(config).get('remap') or {}
    if r.get('bwa_index'):
        return r['bwa_index']
    if starts_at_wasp(config):
        raise ValueError(
            "'workflow: start: wasp' of ATAC bams needs 'modules: wasp: atac: "
            "remap: bwa_index' (the bwa-indexed fasta of the genome the bams "
            "were aligned to, e.g. <cellranger-arc reference>/fasta/genome.fa)")
    aligner = str((config.get('workflow') or {}).get('aligner')).lower()
    if aligner != 'cellranger_arc':
        raise ValueError(
            "WASP of ATAC bams needs 'modules: wasp: atac: remap: bwa_index', "
            "or a multiome run ('workflow: aligner: cellranger_arc')")
    genome = config['modules']['cellranger_arc']['ref']['genome']
    return f"{genome.rstrip('/')}/fasta/genome.fa"


def bwa_extra(config) -> str:
    """bwa mem options ('atac: remap: extra_params'; default -M -Y, as the
    lab's runs)."""
    r = atac_config(config).get('remap') or {}
    extra = r.get('extra_params')
    return "-M -Y" if extra is None else str(extra)


def drop_duplicates(config) -> bool:
    """Whether the duplicates marked by cellranger-arc are removed."""
    mode = str(atac_config(config).get('drop_duplicates') or 'marked').lower()
    if mode not in DUPLICATE_MODES:
        raise ValueError(
            f"'modules: wasp: atac: drop_duplicates' is '{mode}': expected one "
            f"of {', '.join(DUPLICATE_MODES)}")
    return mode == 'marked'
