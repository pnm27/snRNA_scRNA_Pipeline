"""Path of the per-donor mito read counts (split_bams, 'gt_check: no')."""

import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from lib.config_defaults import fill_defaults  # noqa: E402
from lib.params import get_mito_file  # noqa: E402

WC = SimpleNamespace(pool='samp1', donor='D1')


def split_config(**split):
    config = {'workflow': {'start': 'split_bams'},
              'modules': {'split_bams': split}}
    return fill_defaults(config)


def test_next_to_per_donor_bams():
    config = split_config(split_bams_dir='out/', fold_struct_bams='Sample_{pool}-cDNA/')
    assert get_mito_file(WC, config) == 'out/Sample_samp1-cDNA/samp1_mito_reads_stats.txt'


def test_none_with_gt_check_or_null():
    config = split_config()
    config['workflow']['gt_check'] = True
    assert get_mito_file(WC, config) is False
    config = split_config(mito_reads_file=None)
    assert get_mito_file(WC, config) is False
