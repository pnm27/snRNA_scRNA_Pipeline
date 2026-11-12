"""Bams are split by the cDNA demultiplexing: '{modality}' of the split_bams
folder structures is filled with cDNA (lib/utils.py: apply_split_modality),
and a multiome h5ad is split by its cDNA vireo column."""

from lib.utils import apply_split_modality
from lib.params import _split_h5ad_column


def test_modality_filled_in_split_folds():
    config = {'modules': {'split_bams': {
        'fold_struct_hash': 'Sample-{pool}-{modality}',
        'fold_struct_bams': 'Sample_{pool}-{modality}/',
        'fold_struct_donor': '{donor}'}}}
    apply_split_modality(config)
    sb = config['modules']['split_bams']
    assert sb['fold_struct_hash'] == 'Sample-{pool}-cDNA'
    assert sb['fold_struct_bams'] == 'Sample_{pool}-cDNA/'
    assert sb['fold_struct_donor'] == '{donor}'


def test_folds_without_modality_unchanged():
    config = {'modules': {'split_bams': {
        'fold_struct_hash': 'Sample-{pool}-cDNA', 'fold_struct_bams': 'Sample_{pool}-cDNA/'}}}
    apply_split_modality(config)
    assert config['modules']['split_bams']['fold_struct_hash'] == 'Sample-{pool}-cDNA'


def test_h5ad_split_column():
    split_by = {'input': 'h5ad', 'demux': 'vireo', 'column': None}
    assert _split_h5ad_column(split_by) == 'SubID_vs'
    assert _split_h5ad_column(split_by, multiome=True) == 'SubID_vs_cdna'
    # A column set in the config, and the final h5ad's, are kept
    assert _split_h5ad_column({**split_by, 'column': 'x'}, multiome=True) == 'x'
    assert _split_h5ad_column({'input': 'final'}, multiome=True) == 'Donor_final'


def test_mito_file_has_no_placeholders():
    # a params function's path isn't formatted with the wildcards, so the
    # '{modality}' of fold_struct_bams must be filled there too
    from types import SimpleNamespace
    from lib.params import get_mito_file

    config = {'workflow': {'gt_check': False}, 'modules': {'split_bams': {
        'split_bams_dir': 'split/final/', 'fold_struct_bams': 'Sample_{pool}-{modality}/',
        'mito_reads_file': '_mito_reads_stats.txt'}}}
    path = get_mito_file(SimpleNamespace(pool='P1', donor='D1'), config)
    assert path == 'split/final/Sample_P1-cDNA/P1_mito_reads_stats.txt'
