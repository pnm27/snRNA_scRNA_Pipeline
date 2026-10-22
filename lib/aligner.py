#!/usr/bin/env python3
"""Outputs of the aligner of this run.

The prefix of 'last_step' picks the aligner:
- starsolo_*   : STARsolo (cDNA)
- cellranger_* : cellranger count (cDNA)
- multiome_*   : cellranger-arc count (cDNA and ATAC)

Rules and input functions ask this module for the count matrix, bams and
other alignment outputs instead of building the paths of one aligner.
All paths are templates with {pool} (and the wildcards of 'fold_struct').
"""

from types import SimpleNamespace
from .modality import CELLRANGER_ARC


# Keys of the bam and bai of STARsolo and cellranger count (cDNA only)
BAM_KEYS = SimpleNamespace(bam='bam', bai='bai')


# last_step prefix -> aligner
ALIGNERS = {
    'starsolo': 'starsolo',
    'cellranger': 'cellranger',
    'multiome': 'cellranger_arc',
}

# Config section of each aligner
SECTIONS = {
    'starsolo': 'STARsolo_pipeline',
    'cellranger': 'cellranger_count',
    'cellranger_arc': 'cellranger_arc_count',
}


def get_aligner(config) -> str:
    step = config['last_step'].lower()
    for prefix, aligner in ALIGNERS.items():
        if step.startswith(prefix):
            return aligner
    # e.g. 'all_multi_vcf'
    return 'starsolo'


def out_prefix(config, aligner: str | None = None) -> str:
    """Prefix of the outputs of a pool (aligner: this run's, by default)."""
    aligner = aligner or get_aligner(config)
    if aligner == 'starsolo':
        return f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}"
    return f"{config[SECTIONS[aligner]]['bams_dir']}{{pool}}/"


def barcodes_have_suffix(config) -> bool:
    """Whether cell barcodes end with '-1' (cellranger)."""
    return get_aligner(config) != 'starsolo'


# STARsolo's count matrix ('STARsolo_pipeline: count_matrix') -> its
# matrix, barcodes and features keys in 'STARsolo_pipeline'
STARSOLO_MATRIX_KEYS = {
    'GeneFull': ('genefull_matrix', 'genefull_barcodes', 'genefull_features'),
    'Gene': ('gene_matrix', 'gene_barcodes', 'gene_features'),
}


def starsolo_matrix_keys(config) -> tuple:
    return STARSOLO_MATRIX_KEYS[config['STARsolo_pipeline'].get('count_matrix', 'GeneFull')]


def count_matrix(config) -> str:
    """Filtered count matrix (matrix.mtx.gz) of a pool."""
    if get_aligner(config) != 'starsolo':
        return f"{out_prefix(config)}filtered_feature_bc_matrix/matrix.mtx.gz"

    c = config['STARsolo_pipeline']
    key = starsolo_matrix_keys(config)[0]
    # The scripts read the dir of 'matrix.mtx.gz' (scanpy's read_10x_mtx)
    if not c[key].endswith('matrix.mtx.gz'):
        raise ValueError(
            f"'STARsolo_pipeline: {key}' should end with 'matrix.mtx.gz' "
            f"(it's gzipped after STARsolo runs), not: {c[key]}"
        )
    return f"{out_prefix(config)}{c[key]}"


def bam_paths(config, modality: str, suffix: str | None = None,
    aligner: str | None = None) -> dict:
    """bam and bai of 'modality'.

    suffix: a derived bam, e.g. config['split_bams_pipeline']['short_bam']
    (subset to a chromosome) or ['filt_bam'] (filtered by cell barcode).
    STARsolo's are named <prefix><suffix>; cellranger's replace '.bam' of
    the bam's name with the suffix.
    aligner: this run's, by default.
    """
    aligner = aligner or get_aligner(config)
    prefix = out_prefix(config, aligner)
    c = config[SECTIONS[aligner]]
    # Keys of the bam and bai in the aligner's config section
    keys = CELLRANGER_ARC[modality] if aligner == 'cellranger_arc' else BAM_KEYS
    if suffix is None:
        return {'bam': f"{prefix}{c[keys.bam]}", 'bai': f"{prefix}{c[keys.bai]}"}
    if aligner == 'starsolo':
        bam = f"{prefix}{suffix}"
    else:
        bam = f"{prefix}{c[keys.bam]}".replace('.bam', suffix)
    return {'bam': bam, 'bai': f"{bam}.bai"}
