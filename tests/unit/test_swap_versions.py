"""Demultiplex versions of the swap correction file ('swap_correction:
versions', 'status_column') and 'start: split_bams' with 'input: final'."""

import copy, json, os, sys
from types import SimpleNamespace
import pytest
import yaml

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, REPO)
from lib.config_check import config_problems  # noqa: E402
from lib.config_defaults import fill_defaults  # noqa: E402
from lib.config_models import final_sources, swap_pools, is_kept  # noqa: E402
from lib.params import get_params_demux_final  # noqa: E402
from lib.utils import apply_run_tag  # noqa: E402

SCHEMA = yaml.safe_load(open(os.path.join(REPO, "schema", "config.schema.yaml")))

SWAP = (
    "Sample\tDonor_name\tSolution\tperc_het_consistency\tfinal_status\tVersion_used\n"
    "Sample1\tdonor2\tdonA\t0.95\tTrue\tw_gt\n"
    "Sample1\tdonor3\tdonB\t0.95\tFalse\tw_gt\n"
    "Sample2\tdonor4\tdonC\t0.96\tTrue\two_gt\n"
    "Sample3\tdonor5\tdonD\t0.96\tFalse\tw_gt\n"
)


def make_config(tmp_path, versions, status='final_status', tagged=True):
    swap = tmp_path / "swap.tsv"
    swap.write_text(SWAP)
    config = {
        'workflow': {'start': 'split_bams', 'until': 'split_bams'},
        'modules': {
            'demultiplex': {'swap_correction': {
                'file': str(swap), 'pool_column': 'Sample',
                'version_column': 'Version_used', 'donor_from': 'Donor_name',
                'donor_to': 'Solution', 'status_column': status,
                'versions': versions}},
            'split_bams': {'split_by': {'demux': 'vireo', 'input': 'final'}},
        },
    }
    fill_defaults(config)
    if tagged:  # as the Snakefile does after checking the config
        apply_run_tag(config)
    return config


def test_minimal_final_config_is_valid(tmp_path):
    config = make_config(tmp_path, {'w_gt': 'vireo:w_gt', 'wo_gt': 'vireo:wo_gt'},
                         tagged=False)
    assert config_problems(config, SCHEMA)[0] == []


def test_final_start_needs_swap_file():
    config = fill_defaults({'workflow': {'start': 'split_bams', 'until': 'split_bams'},
                            'modules': {'split_bams': {'split_by': {'input': 'final'}}}})
    errors = config_problems(config, SCHEMA)[0]
    assert any("swap_correction: file" in e for e in errors)


def test_versions_from_json(tmp_path):
    template = tmp_path / "template.json"
    template.write_text(json.dumps({'demultiplex_paths': {
        'w_gt': {'demux': 'runs/w_gt/', 'info': 'info/w_gt'},
        'wo_gt': {'demux': 'runs/wo_gt/'}}}))
    config = make_config(tmp_path, str(template))
    [run] = final_sources(config, 'Sample1')
    assert (run.tool, run.label) == ('vireo', 'vireo:w_gt')
    assert run.h5ad == 'runs/w_gt/Sample_{pool}-cDNA_STARsolo_vS_out.h5ad'
    assert run.stats == 'info/w_gt/Sample_{pool}-cDNA_vireo_stats.tsv'
    [run] = final_sources(config, 'Sample2')
    assert run.stats is None  # no 'info'


def test_versions_map(tmp_path):
    config = make_config(tmp_path, {
        'w_gt': 'both:w_gt',
        'wo_gt': {'demux': 'elsewhere', 'tool': 'solo'}})
    assert [r.label for r in final_sources(config, 'Sample1')] == ['solo:w_gt', 'vireo:w_gt']
    [run] = final_sources(config, 'Sample2')
    assert run.tool == 'solo' and run.h5ad.startswith('elsewhere/')


def test_unknown_version(tmp_path):
    config = make_config(tmp_path, {'w_gt': 'vireo:w_gt'})
    with pytest.raises(ValueError, match="'wo_gt' .* isn't in 'swap_correction: versions'"):
        final_sources(config, 'Sample2')


def test_without_versions_needs_tool_labels(tmp_path):
    config = make_config(tmp_path, None)
    with pytest.raises(ValueError, match="Unknown tool 'w_gt'"):
        final_sources(config, 'Sample1')


def test_pools_with_a_donor_kept(tmp_path):
    config = make_config(tmp_path, None)
    assert swap_pools(config) == ['Sample1', 'Sample2']
    assert swap_pools(config, kept_only=False) == ['Sample1', 'Sample2', 'Sample3']
    config['modules']['demultiplex']['swap_correction']['status_column'] = None
    assert swap_pools(config) == ['Sample1', 'Sample2', 'Sample3']


def test_missing_column(tmp_path):
    config = make_config(tmp_path, None, status='status')
    with pytest.raises(ValueError, match="no column.* status"):
        swap_pools(config)


@pytest.mark.parametrize("value, kept", [
    (True, True), ("True", True), ("yes", True), (1, True), ("1", True),
    (False, False), ("False", False), ("", False), (float('nan'), False)])
def test_is_kept(value, kept):
    assert is_kept(value) is kept


def test_demux_final_gets_status_column(tmp_path):
    config = make_config(tmp_path, {'w_gt': 'vireo:w_gt', 'wo_gt': 'vireo:wo_gt'})
    params = get_params_demux_final(SimpleNamespace(pool='Sample1'),
                                    SimpleNamespace(swap='swap.tsv'), config)
    assert "--swap_status final_status" in params
    assert "--sources 'vireo:w_gt'" in params


def test_demux_final_drops_donors_not_true(tmp_path):
    pytest.importorskip("anndata")
    sys.path.insert(0, os.path.join(REPO, "helper_py_scripts"))
    from demux_final import swap_donor_map
    swap = tmp_path / "swap.tsv"
    swap.write_text(SWAP)
    cols = ['Sample', 'Donor_name', 'Solution']
    m = swap_donor_map(str(swap), cols, 'Sample1', 'final_status')
    assert m['donor2'] == 'donA' and m['donor3'] != m['donor3']  # NaN: dropped
    assert swap_donor_map(str(swap), cols, 'Sample1')['donor3'] == 'donB'
