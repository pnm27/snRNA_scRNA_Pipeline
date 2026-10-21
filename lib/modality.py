#!/usr/bin/env python3
"""Alignment outputs of each modality ('modality' wildcard).

Multiome runs ('last_step' starts with 'multiome') are aligned by
cellranger-arc, which writes one bam per modality (cDNA and ATAC) for each
pool. Other runs are aligned by STARsolo (cDNA only).

The bams derived from these (subset to a chromosome, or filtered by cell
barcode) are named by a suffix from 'split_bams_pipeline' (see bam_paths).
"""

from types import SimpleNamespace


# Per modality of a cellranger-arc run: keys of its bam and bai in
# config['cellranger_arc_count'], and whether its reads carry UMIs
CELLRANGER_ARC = {
    'cDNA': SimpleNamespace(bam='gex_bam', bai='gex_bai', umi=True),
    'ATAC': SimpleNamespace(bam='atac_bam', bai='atac_bai', umi=False),
}


def is_multiome(config) -> bool:
    return 'multiome' in config['last_step'].lower()


def bam_modalities(config) -> list[str]:
    """Modalities with a bam of their own."""
    return list(CELLRANGER_ARC) if is_multiome(config) else ['cDNA']


def get_modality(wildcards) -> str:
    """'modality' wildcard (cDNA, if the rule has none)."""
    modality = getattr(wildcards, 'modality', None) or 'cDNA'
    for m in CELLRANGER_ARC:
        if m.lower() == modality.lower():
            return m
    raise ValueError(f"No bam for the modality '{modality}'")


def bam_paths(config, modality: str, suffix: str | None = None) -> dict:
    """bam and bai of 'modality' (path templates with {pool}).

    suffix: a derived bam, e.g. config['split_bams_pipeline']['short_bam']
    (subset to a chromosome) or ['filt_bam'] (filtered by cell barcode).
    """
    if is_multiome(config):
        c = config['cellranger_arc_count']
        keys = CELLRANGER_ARC[modality]
        bam = f"{c['bams_dir']}{{pool}}/{c[keys.bam]}"
        if suffix is None:
            return {'bam': bam, 'bai': f"{c['bams_dir']}{{pool}}/{c[keys.bai]}"}
        bam = bam.replace('.bam', suffix)
    else:
        c = config['STARsolo_pipeline']
        prefix = f"{c['bams_dir']}{config['fold_struct']}"
        if suffix is None:
            return {'bam': f"{prefix}{c['bam']}", 'bai': f"{prefix}{c['bai']}"}
        bam = f"{prefix}{suffix}"
    return {'bam': bam, 'bai': f"{bam}.bai"}


def has_umi(config, modality: str) -> bool:
    return CELLRANGER_ARC[modality].umi if is_multiome(config) else True
