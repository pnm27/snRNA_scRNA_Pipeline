#!/usr/bin/env python3
"""Outputs of the aligner of this run.

'workflow: aligner' picks the aligner:
- STARsolo       : STARsolo (cDNA)
- cellranger     : cellranger count (cDNA)
- cellranger_arc : cellranger-arc count (cDNA and ATAC)

Rules and input functions ask this module for the count matrix, bams and
other alignment outputs instead of building the paths of one aligner.
All paths are templates with {pool} (and the wildcards of 'fold_struct').
"""

from types import SimpleNamespace
from .modality import CELLRANGER_ARC
from .utils import run_settings


# Keys of the bam and bai of STARsolo and cellranger count (cDNA only)
BAM_KEYS = SimpleNamespace(bam='bam', bai='bai')


# Module (config['modules']) of each aligner
SECTIONS = {
    'starsolo': 'STARsolo',
    'cellranger': 'cellranger',
    'cellranger_arc': 'cellranger_arc',
}


def get_aligner(config) -> str:
    """'workflow: aligner', lowercase (a key of SECTIONS)."""
    aligner = run_settings(config).aligner
    if aligner not in SECTIONS:
        raise ValueError(f"Unrecognised aligner: '{config['workflow']['aligner']}'")
    return aligner


def out_prefix(config, aligner: str | None = None) -> str:
    """Prefix of the outputs of a pool (aligner: this run's, by default)."""
    aligner = aligner or get_aligner(config)
    if aligner == 'starsolo':
        return f"{config['modules']['STARsolo']['bams_dir']}{config['modules']['STARsolo']['fold_struct']}"
    return f"{config['modules'][SECTIONS[aligner]]['bams_dir']}{{pool}}/"


def barcodes_have_suffix(config) -> bool:
    """Whether cell barcodes end with '-1' (cellranger)."""
    return get_aligner(config) != 'starsolo'


# STARsolo's count matrix ('modules: STARsolo: count_matrix') -> its
# matrix, barcodes and features keys in 'modules: STARsolo'
STARSOLO_MATRIX_KEYS = {
    'GeneFull': ('genefull_matrix', 'genefull_barcodes', 'genefull_features'),
    'Gene': ('gene_matrix', 'gene_barcodes', 'gene_features'),
}


def starsolo_matrix_keys(config) -> tuple:
    return STARSOLO_MATRIX_KEYS[config['modules']['STARsolo'].get('count_matrix', 'GeneFull')]


def count_matrix(config) -> str:
    """Filtered count matrix (matrix.mtx.gz) of a pool."""
    if get_aligner(config) != 'starsolo':
        return f"{out_prefix(config)}filtered_feature_bc_matrix/matrix.mtx.gz"

    c = config['modules']['STARsolo']
    key = starsolo_matrix_keys(config)[0]
    # The scripts read the dir of 'matrix.mtx.gz' (scanpy's read_10x_mtx)
    if not c[key].endswith('matrix.mtx.gz'):
        raise ValueError(
            f"'modules: STARsolo: {key}' should end with 'matrix.mtx.gz' "
            f"(it's gzipped after STARsolo runs), not: {c[key]}"
        )
    return f"{out_prefix(config)}{c[key]}"


def bam_paths(config, modality: str, suffix: str | None = None,
    aligner: str | None = None) -> dict:
    """bam and bai of 'modality'.

    suffix: a derived bam, e.g. config['modules']['split_bams']['short_bam']
    (subset to a chromosome) or ['filt_bam'] (filtered by cell barcode).
    STARsolo's are named <prefix><suffix>; cellranger's replace '.bam' of
    the bam's name with the suffix.
    aligner: this run's, by default.
    """
    aligner = aligner or get_aligner(config)
    prefix = out_prefix(config, aligner)
    c = config['modules'][SECTIONS[aligner]]
    # Keys of the bam and bai in the aligner's config section
    keys = CELLRANGER_ARC[modality] if aligner == 'cellranger_arc' else BAM_KEYS
    if suffix is None:
        return {'bam': f"{prefix}{c[keys.bam]}", 'bai': f"{prefix}{c[keys.bai]}"}
    if aligner == 'starsolo':
        bam = f"{prefix}{suffix}"
    else:
        bam = f"{prefix}{c[keys.bam]}".replace('.bam', suffix)
    return {'bam': bam, 'bai': f"{bam}.bai"}
