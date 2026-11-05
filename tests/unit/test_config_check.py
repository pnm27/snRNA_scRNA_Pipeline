"""lib/config_check.py: each schema rule catches its mistake and says why
(needs jsonschema; run with pytest)."""

import copy, os, sys
import pytest
import yaml

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, REPO)
from lib.config_check import ConfigError, check_config, config_problems, duplicate_keys

SCHEMA_PATH = os.path.join(REPO, "schema", "config.schema.yaml")
SCHEMA = yaml.safe_load(open(SCHEMA_PATH))
EXAMPLE = yaml.safe_load(open(os.path.join(REPO, "new_config.example.yaml")))


def merge(base, over):
    out = copy.deepcopy(base)
    for k, v in over.items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def demux(**kw):
    return {'experiment': {'demultiplexing': kw}}


def external(**kw):
    return {'modules': {'split_bams': {'split_by': {'input': 'external', 'external': kw}}}}


START_SPLIT = {'workflow': {'start': 'split_bams', 'until': 'split_bams'}}


def test_example_is_valid():
    assert config_problems(EXAMPLE, SCHEMA) == ([], [])


# override of the example -> (text of the error, text of its explanation)
CASES = {
    'unknown key': (
        {'workflow': {'untill': 'demux'}},
        "Unknown key(s) 'untill' in workflow", "migrate_config.py"),
    'old key at the top level': (
        {'run_tag': 'w_gt'},
        "Unknown key(s) 'run_tag' in (top level)", "only these keys are read"),
    'until typo': (
        {'workflow': {'until': 'demultiplex'}},
        "workflow: until: 'demultiplex' is not one of", "Last stage to run"),
    'method typo': (
        demux(method='hashsolo'),
        "method: 'hashsolo' is not one of", "Tool run now"),
    'add to solo with solo': (
        demux(method='solo', add_to_previous='solo'),
        "method: 'vireo' was expected", "means method: vireo"),
    'add to vireo with vireo': (
        merge(demux(method='vireo', add_to_previous='vireo'),
              {'modules': {'demultiplex': {'previous_run_tag': 'w_gt'}}}),
        "method: 'solo' was expected", "means method: solo"),
    'add without previous_run_tag': (
        demux(method='vireo', add_to_previous='solo'),
        "previous_run_tag: None is not of type 'string'", "found by its run_tag"),
    'multi-VCF with hashsolo': (
        merge(demux(method='both', genotype={'enabled': True, 'multi_vcf': True}),
              {'modules': {'gt_demux': {'fold_struct': 'Sample_{pool}-cDNA/{vcf_type}/'}}}),
        "method: 'vireo' was expected", "support only vireo"),
    'multi-VCF without the vcf_type wildcard': (
        demux(method='vireo', genotype={'enabled': True, 'multi_vcf': True}),
        "fold_struct: 'Sample_{pool}-cDNA/' does not match", "'vcf_type' wildcard must be in"),
    'multi-VCF split bams': (
        merge(demux(method='vireo', genotype={'enabled': True, 'multi_vcf': True}),
              {'workflow': {'until': 'split_bams'},
               'modules': {'gt_demux': {'fold_struct': 'Sample_{pool}-cDNA/{vcf_type}/'}}}),
        "until: 'split_bams' isn't allowed here", "isn't supported for multi-VCF runs"),
    'split bams without demultiplexing': (
        merge(demux(method=None), {'workflow': {'until': 'split_bams'}}),
        "method: None is not one of ['solo', 'vireo', 'both']", "needs a demultiplexing method"),
    'hashsolo with multiome': (
        {'workflow': {'aligner': 'cellranger_arc'}},
        "method: 'solo' isn't allowed here", "isn't supported for multiome runs yet"),
    'hashsolo with cellranger count': (
        {'workflow': {'aligner': 'cellranger', 'picard_metrics': None}},
        "method: 'solo' isn't allowed here", "isn't supported with cellranger count yet"),
    'cellranger count without its module': (
        merge(demux(method='vireo'),
              {'workflow': {'aligner': 'cellranger', 'picard_metrics': None},
               'modules': {'cellranger': None}}),
        "modules: cellranger: None is not of type 'object'", "need 'modules: cellranger'"),
    'hashsolo without HTO': (
        {'experiment': {'hto': {'enabled': False}}},
        "hto: enabled: True was expected", "demultiplexes by the HTO counts"),
    'hashsolo without the wet lab file': (
        {'experiment': {'inputs': {'wet_lab_info': None}}},
        "wet_lab_info: None is not of type 'string'", "maps hashtags to samples"),
    'no start': (
        {'workflow': {'start': None}},
        "start: None is not of type 'string'", "First stage to run"),
    'start split_bams without external labels': (
        START_SPLIT,
        "input: 'raw' does not match", "splits bams that are already there"),
    'start split_bams until demux': (
        merge({'workflow': {'start': 'split_bams', 'until': 'demux'}},
              external(file='p.h5ad', pool_column='pool')),
        "until: 'demux' is not one of", "splits bams that are already there"),
    'external without file': (
        external(pool_column='pool'),
        "external: file: None is not of type 'string'", "needs 'split_by: external: file'"),
    'external without pool_column': (
        external(file='p.h5ad'),
        "pool_column: None is not of type 'string'", "needs 'split_by: external: file'"),
    'external key typo': (
        external(file='p.h5ad', pool_column='pool', pool_col='x'),
        "Unknown key(s) 'pool_col' in modules: split_bams: split_by: external", "only these keys"),
    'external step without its argument': (
        external(file='p.h5ad', pool_column='pool',
                 barcode_transformations=[{'type': 'pythonic', 'operation': 'split'}]),
        "is not valid under any of the given schemas", "Each step is"),
    'pools table with start align': (
        external(file='p.h5ad', pool_column='pool', pools='bams.tsv'),
        "start: 'split_bams' was expected", "lists bams that are already there"),
    'pools table on a multiome run': (
        merge(merge(START_SPLIT, {'workflow': {'aligner': 'cellranger_arc'}}),
              external(file='p.h5ad', pool_column='pool', pools='bams.tsv')),
        "aligner", "isn't supported for multiome runs"),
    'HTO without multiplexing': (
        merge(demux(method='vireo'), {'experiment': {'hto': {'enabled': True, 'multiplexing': None}}}),
        "multiplexing: None is not one of ['single', 'multiple']", "picks the kite rules"),
}


@pytest.mark.parametrize("override, error, why", CASES.values(), ids=CASES.keys())
def test_rule(override, error, why):
    errors, _ = config_problems(merge(EXAMPLE, override), SCHEMA)
    assert any(error in e and why in e for e in errors), errors


def test_picard_warning():
    config = merge(EXAMPLE, merge(demux(method='vireo'),
                                  {'workflow': {'aligner': 'cellranger', 'picard_metrics': 'all'}}))
    errors, warnings = config_problems(config, SCHEMA)
    assert errors == [] and any("Picard runs on STARsolo runs only" in w for w in warnings)


def test_all_errors_at_once():
    config = merge(EXAMPLE, {'workflow': {'untill': 'x', 'until': 'y'}})
    with pytest.raises(ConfigError, match=r"2 problem\(s\)"):
        check_config(config, SCHEMA_PATH)


def test_old_layout():
    with pytest.raises(ConfigError, match="old layout"):
        check_config({'last_step': 'STARsolo_demux'}, SCHEMA_PATH)


def test_duplicate_keys(tmp_path):
    f = tmp_path / "c.yaml"
    f.write_text("a: &x\n  k: 1\nb:\n  <<: *x\n  j: 1\n  j: 2\n")
    msg = duplicate_keys(str(f))
    assert "Duplicate key 'j' (lines 5 and 6)" in msg and "YAML keeps the last value" in msg
    f.write_text("a: &x\n  k: 1\nb:\n  <<: *x\n  k: 2\n")  # overriding a merged key is fine
    assert duplicate_keys(str(f)) is None


def test_external_needs_no_demux():
    """Bams split by external labels need no demultiplexing method."""
    config = merge(EXAMPLE, merge(merge(START_SPLIT, demux(method=None)),
                                  external(file='p.h5ad', pool_column='pool')))
    assert config_problems(config, SCHEMA) == ([], [])


def test_step_defs_match_h5ad_schema():
    """The steps of split_by: external are those of the h5ad annotations."""
    import json
    with open(os.path.join(REPO, "schema", "demul_samples_annotate_h5ad.schema.json")) as f:
        assert SCHEMA['$defs'] == json.load(f)['$defs']


def test_pools_table_with_start_split_bams():
    """A pools table is fine for bams that are already there."""
    config = merge(EXAMPLE, merge(merge(START_SPLIT, demux(method=None)),
                                  external(file='p.h5ad', pool_column='pool',
                                           pools='bams.tsv')))
    assert config_problems(config, SCHEMA) == ([], [])
