#!/usr/bin/env python3
"""Modalities ('modality' wildcard) of a run.

Multiome runs are aligned by cellranger-arc (workflow: aligner:
cellranger_arc), which writes one bam per modality (cDNA and ATAC) for each
pool. Other runs have only cDNA. Their bams are in lib/aligner.py.
"""

from types import SimpleNamespace
from .utils import run_settings, SPLIT_MODALITY


# Per modality of a cellranger-arc run: keys of its bam and bai in
# config['modules']['cellranger_arc'], and whether its reads carry UMIs
CELLRANGER_ARC = {
    'cDNA': SimpleNamespace(bam='gex_bam', bai='gex_bai', umi=True),
    'ATAC': SimpleNamespace(bam='atac_bam', bai='atac_bai', umi=False),
}


def is_multiome(config) -> bool:
    return run_settings(config).aligner == 'cellranger_arc'


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


def normalize_modality(name) -> str:
    """'cDNA' or 'ATAC' (any case); else a ValueError."""
    for m in CELLRANGER_ARC:
        if m.lower() == str(name).strip().lower():
            return m
    raise ValueError(f"Unknown modality '{name}' (expected one of: {', '.join(CELLRANGER_ARC)})")


def split_modalities(config) -> list[str]:
    """Modalities whose bams are split per donor ('modules: split_bams:
    modalities', default cDNA). Every modality is split by the same donor
    labels: the cDNA demultiplexing, or the external labels ('split_by:
    input: external'). ATAC needs a multiome run (aligner cellranger_arc)."""
    sb = (config.get('modules') or {}).get('split_bams') or {}
    given = sb.get('modalities') or [SPLIT_MODALITY]
    if isinstance(given, str):
        given = [given]
    out = []
    for m in given:
        m = normalize_modality(m)
        if m not in out:
            out.append(m)
    if 'ATAC' in out and not is_multiome(config):
        raise ValueError(
            "'modules: split_bams: modalities' has ATAC: only multiome runs "
            "('workflow: aligner: cellranger_arc') have an ATAC bam")
    return out
