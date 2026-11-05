"""lib/config_defaults.py: a config of existing bams ('start: split_bams')
needs only what differs from the defaults (needs jsonschema)."""

import copy, os, sys
import yaml

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, REPO)
from lib.config_check import config_problems
from lib.config_defaults import fill_defaults, merge

SCHEMA = yaml.safe_load(open(os.path.join(REPO, "schema", "config.schema.yaml")))

MINIMAL = {
    'workflow': {'start': 'split_bams', 'until': 'split_bams'},
    'modules': {'split_bams': {'split_by': {'input': 'external', 'external': {
        'file': 'cellBC_vs_origin.txt', 'barcode_column': 'barcode',
        'pool_column': 'pool', 'pools': 'bams.tsv'}}}},
}


def test_minimal_config_is_valid():
    config = fill_defaults(copy.deepcopy(MINIMAL))
    assert config_problems(config, SCHEMA) == ([], [])
    assert config['workflow']['gt_check'] is False
    assert config['experiment']['demultiplexing']['method'] is None


def test_config_wins_over_defaults():
    over = merge(MINIMAL, {
        'workflow': {'gt_check': True},
        'modules': {'split_bams': {
            'split_bams_dir': 'out/',
            'split_by': {'external': {'barcode_transformations': []}}}}})
    config = fill_defaults(copy.deepcopy(over))
    assert config['workflow']['gt_check'] is True
    assert config['modules']['split_bams']['split_bams_dir'] == 'out/'
    ext = config['modules']['split_bams']['split_by']['external']
    assert ext['barcode_transformations'] == [] and ext['donor_column'] == 'individualID'


def test_defaults_keep_errors():
    """A mistake in the config is still reported."""
    config = fill_defaults(merge(MINIMAL, {'workflow': {'untill': 'x'}}))
    errors, _ = config_problems(config, SCHEMA)
    assert any("Unknown key(s) 'untill'" in e for e in errors)


def test_align_config_unchanged():
    config = {'workflow': {'start': 'align', 'until': 'align'}}
    assert fill_defaults(config) == {'workflow': {'start': 'align', 'until': 'align'}}
    assert config_problems(config, SCHEMA)[0]  # still missing its sections
