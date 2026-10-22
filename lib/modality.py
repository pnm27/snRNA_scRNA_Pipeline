#!/usr/bin/env python3
"""Modalities ('modality' wildcard) of a run.

Multiome runs are aligned by cellranger-arc (workflow: aligner:
cellranger_arc), which writes one bam per modality (cDNA and ATAC) for each
pool. Other runs have only cDNA. Their bams are in lib/aligner.py.
"""

from types import SimpleNamespace


# Per modality of a cellranger-arc run: keys of its bam and bai in
# config['modules']['cellranger_arc'], and whether its reads carry UMIs
CELLRANGER_ARC = {
    'cDNA': SimpleNamespace(bam='gex_bam', bai='gex_bai', umi=True),
    'ATAC': SimpleNamespace(bam='atac_bam', bai='atac_bai', umi=False),
}


def is_multiome(config) -> bool:
    return config['workflow']['aligner'].lower() == 'cellranger_arc'


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


def has_umi(config, modality: str) -> bool:
    return CELLRANGER_ARC[modality].umi if is_multiome(config) else True
