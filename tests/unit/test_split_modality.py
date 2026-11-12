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
