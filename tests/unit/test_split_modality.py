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


# ATAC bams of multiome, split by the same donor labels ('split_bams: modalities')

import pytest
from lib.utils import split_bams_fold
from lib.modality import split_modalities
from lib.io import barcode_check_path, external_barcode_check


def _config(fold, modalities=None, aligner='cellranger_arc', input='raw'):
    c = {'workflow': {'aligner': aligner, 'gt_check': False},
         'modules': {'split_bams': {
             'fold_struct_hash': 'Sample-{pool}-{modality}', 'fold_struct_bams': fold,
             'fold_struct_donor': '{donor}', 'inp_split_bams_dir': 'inp/',
             'split_bams_dir': 'split/final/', 'mito_reads_file': '_mito.txt',
             'split_by': {'input': input}}}}
    if modalities:
        c['modules']['split_bams']['modalities'] = modalities
    apply_split_modality(c)
    return c


def test_atac_fold():
    # '{modality}' filled with ATAC; without one, an ATAC/ folder in front
    c = _config('Sample_{pool}-{modality}/')
    assert split_bams_fold(c) == 'Sample_{pool}-cDNA/'
    assert split_bams_fold(c, 'ATAC') == 'Sample_{pool}-ATAC/'
    c = _config('Sample_{pool}-cDNA/')
    assert split_bams_fold(c, 'ATAC') == 'ATAC/Sample_{pool}-cDNA/'


def test_atac_mito_file():
    from types import SimpleNamespace
    from lib.params import get_mito_file
    c = _config('Sample_{pool}-{modality}/')
    path = get_mito_file(SimpleNamespace(pool='P1', donor='D1'), c, modality='ATAC')
    assert path == 'split/final/Sample_P1-ATAC/P1_mito.txt'


def test_split_modalities():
    assert split_modalities(_config('x/')) == ['cDNA']
    assert split_modalities(_config('x/', ['cdna', 'ATAC', 'atac'])) == ['cDNA', 'ATAC']
    with pytest.raises(ValueError, match="only multiome runs"):
        split_modalities(_config('x/', ['ATAC'], aligner='STARsolo'))
    with pytest.raises(ValueError, match="Unknown modality"):
        split_modalities(_config('x/', ['HTO']))


def test_atac_barcode_check():
    # the ATAC bam's CB tags are always checked against the (cDNA) labels
    c = _config('x/', ['cDNA', 'ATAC'])
    assert external_barcode_check(config=c) == []
    assert external_barcode_check(config=c, modality='ATAC') == [
        'inp/Sample-{pool}-cDNA_bc_check_ATAC.txt']
    c = _config('x/', input='external')
    assert external_barcode_check(config=c) == [barcode_check_path(c)]
