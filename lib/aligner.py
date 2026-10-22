#!/usr/bin/env python3
"""Outputs of the aligner of this run.

The prefix of 'last_step' picks the aligner:
- starsolo_*  : STARsolo (cDNA)
- multiome_*  : cellranger-arc count (cDNA and ATAC)

Rules and input functions ask this module for the count matrix, bams and
other alignment outputs instead of building the paths of one aligner.
All paths are templates with {pool} (and the wildcards of 'fold_struct').
"""

from .modality import CELLRANGER_ARC


# last_step prefix -> aligner
ALIGNERS = {
    'starsolo': 'starsolo',
    'multiome': 'cellranger_arc',
}


def get_aligner(config) -> str:
    step = config['last_step'].lower()
    for prefix, aligner in ALIGNERS.items():
        if step.startswith(prefix):
            return aligner
    # e.g. 'all_multi_vcf'
    return 'starsolo'


def out_prefix(config) -> str:
    """Prefix of the aligner's outputs for a pool."""
    if get_aligner(config) == 'starsolo':
        return f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}"
    return f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/"


def barcodes_have_suffix(config) -> bool:
    """Whether cell barcodes end with '-1' (cellranger)."""
    return get_aligner(config) != 'starsolo'


def count_matrix(config) -> str:
    """Filtered count matrix (matrix.mtx.gz) of a pool."""
    if get_aligner(config) == 'starsolo':
        return f"{out_prefix(config)}{config['STARsolo_pipeline']['genefull_matrix']}"
    return f"{out_prefix(config)}filtered_feature_bc_matrix/matrix.mtx.gz"


def bam_paths(config, modality: str, suffix: str | None = None) -> dict:
    """bam and bai of 'modality'.

    suffix: a derived bam, e.g. config['split_bams_pipeline']['short_bam']
    (subset to a chromosome) or ['filt_bam'] (filtered by cell barcode).
    STARsolo's are named <prefix><suffix>; cellranger's replace '.bam' of
    the bam's name with the suffix.
    """
    prefix = out_prefix(config)
    if get_aligner(config) == 'starsolo':
        c = config['STARsolo_pipeline']
        if suffix is None:
            return {'bam': f"{prefix}{c['bam']}", 'bai': f"{prefix}{c['bai']}"}
        bam = f"{prefix}{suffix}"
    else:
        c = config['cellranger_arc_count']
        keys = CELLRANGER_ARC[modality]
        if suffix is None:
            return {'bam': f"{prefix}{c[keys.bam]}", 'bai': f"{prefix}{c[keys.bai]}"}
        bam = f"{prefix}{c[keys.bam]}".replace('.bam', suffix)
    return {'bam': bam, 'bai': f"{bam}.bai"}
