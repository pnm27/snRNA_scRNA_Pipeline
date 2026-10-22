"""migrate_config.py: old -> new config layout (needs ruamel.yaml,
pandas and snakemake; run with pytest).

old_layout/ has the example config and the dry-run scenarios in the old
layout; migrated, they must match the files the pipeline uses now. The
new keys are checked against the old keys' meaning (OLD_* below, as the
pipeline read them before the restructure)."""

import ast, copy, glob, os, sys
import pytest
import yaml

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(REPO, "helper_py_scripts"))
sys.path.insert(0, REPO)
from migrate_config import MOVES, UNUSED, migrate_text
from lib.aligner import get_aligner
from lib.config_models import demux_tools, is_multi_vcf
from lib.utils import is_genotyped

OLD = os.path.join(os.path.dirname(__file__), "old_layout")
EXAMPLE = os.path.join(REPO, "new_config.example.yaml")
OLD_EXAMPLE = os.path.join(OLD, "new_config.example.yaml")
OLD_SCENARIOS = sorted(glob.glob(os.path.join(OLD, "scenarios", "*.yaml")))
# Old keys whose meaning moved to new keys (checked by check_derived)
DERIVED = {'last_step', 'demux_type', 'demux_run_type', 'hto_demux_type'}


# How the pipeline read the old keys ---------------------------------------
OLD_SPECIAL_STEPS = {
    "starsolo_kb_solo":            ["align"],
    "starsolo_cellsnp":            ["cellsnp", "picard"],
    "starsolo_demux":              ["demux", "picard"],
    "starsolo_gt_demux_multi_vcf": ["demux"],
    "all_multi_vcf":               ["align", "demux"],
    "cellranger_cellsnp":          ["cellsnp"],
    "cellranger_demux":            ["demux"],
    "multiome_alignment":          ["align"],
    "multiome_gt_demux":           ["demux"],
}


def old_stages(step):
    step = step.lower()
    if step in OLD_SPECIAL_STEPS:
        return OLD_SPECIAL_STEPS[step]
    stages = ["align"]
    if any(s in step for s in ["gt_demux", "split_bams", "identify_swaps"]):
        stages.append("demux")
    if any(s in step for s in ["split_bams", "identify_swaps"]):
        stages.append("split_bams")
    if "identify_swaps" in step:
        stages.append("identify_swaps")
    return stages


def old_aligner(step):
    step = step.lower()
    for prefix, aligner in [('starsolo', 'starsolo'), ('cellranger', 'cellranger'),
                            ('multiome', 'cellranger_arc')]:
        if step.startswith(prefix):
            return aligner
    return 'starsolo'


OLD_TOOLS = {None: ['none'], 'solo': ['solo'], 'add_solo': ['solo'], 'vireo': ['vireo'],
             'add_vireo': ['vireo'], 'both': ['solo', 'vireo']}


def until_stages(until):
    """rules/produce_targets.smk's UNTIL_STAGES (the .smk can't be imported)."""
    src = open(os.path.join(REPO, "rules", "produce_targets.smk")).read()
    keep = [n for n in ast.parse(src).body
            if isinstance(n, ast.Assign) and n.targets[0].id == 'UNTIL_STAGES']
    ns = {}
    exec(compile(ast.Module(keep, []), "produce_targets.smk", "exec"), ns)
    return ns['UNTIL_STAGES'][until]


def load(path):
    with open(path) as f:
        return f.read()


def migrated(text):
    """The migrated config, as snakemake (PyYAML) reads it."""
    return yaml.safe_load(migrate_text(text))


def pop_path(d, path):
    *parents, key = path.split('.')
    node = d
    for p in parents:
        if p not in node:
            return KeyError
        node = node[p]
    return node.pop(key, KeyError)


def back(new):
    """The old keys moved by MOVES, taken back out of the new layout; the
    derived keys are left in 'new'."""
    old = {}
    # Folder structures first: they are inside the moved blocks
    for key in sorted(MOVES, key=lambda k: not k.startswith('fold_struct')):
        value = pop_path(new, MOVES[key])
        if value is not KeyError:
            old[key] = value
    return old


def merge(base, over):
    """snakemake's merge of --configfile's (nested dicts are merged)."""
    out = copy.deepcopy(base)
    for k, v in over.items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def check_derived(old, new):
    """The new keys mean, to the pipeline, what the old ones did."""
    exp, wf = new.get('experiment', {}), new.get('workflow', {})
    if 'last_step' in old:
        step = old['last_step']
        assert get_aligner(new) == old_aligner(step)
        # 'until' always runs the alignment (Picard's metrics are part of it)
        before = set(old_stages(step)) - {'picard'}
        after = set(until_stages(wf['until']))
        assert before <= after and after - before <= {'align'}
        assert is_multi_vcf(new) == step.lower().endswith('multi_vcf')
    if 'demux_type' in old:
        demux_type = old['demux_type']
        demux_type = None if str(demux_type).lower() in ('none', 'null') else demux_type
        assert demux_tools(new) == OLD_TOOLS[demux_type and demux_type.lower()]
        assert (exp['demultiplexing']['add_to_previous'] is not None) == str(demux_type).startswith('add_')
    if 'demux_run_type' in old:
        assert is_genotyped(new) == (old['demux_run_type'] != 'wo_gt')
    if 'hto_demux_type' in old:
        hto = exp['hto']
        assert hto['enabled'] == (old['hto_demux_type'] is not None)
        if hto['enabled']:
            assert hto['multiplexing'] == old['hto_demux_type'].lower()


def check_round_trip(text):
    old = yaml.safe_load(text)
    new = migrated(text)
    check_derived(old, new)
    expected = {k: v for k, v in old.items() if k not in DERIVED | UNUSED}
    # A block holding only folder structures is empty once they are out
    moved = {k: v for k, v in back(new).items() if v != {} or k in expected}
    assert moved == expected
    # Nothing left but the derived keys
    for path in ['workflow.aligner', 'workflow.until', 'experiment.hto.enabled',
                 'experiment.hto.multiplexing', 'experiment.demultiplexing.method',
                 'experiment.demultiplexing.add_to_previous',
                 'experiment.demultiplexing.genotype.enabled',
                 'experiment.demultiplexing.genotype.multi_vcf']:
        pop_path(new, path)

    def empty(d):
        return all(isinstance(v, dict) and empty(v) for v in d.values())
    assert empty(new), new


def test_example():
    check_round_trip(load(OLD_EXAMPLE))
    assert migrated(load(OLD_EXAMPLE)) == yaml.safe_load(load(EXAMPLE))


@pytest.mark.parametrize("scenario", OLD_SCENARIOS, ids=os.path.basename)
def test_scenario(scenario):
    check_round_trip(load(scenario))
    now = os.path.join(REPO, "tests", "dryrun", "scenarios", os.path.basename(scenario))
    assert migrated(load(scenario)) == yaml.safe_load(load(now))


def test_header_comments():
    text = load(OLD_SCENARIOS[0])
    header = [l for l in text.splitlines() if l.startswith('#')]
    assert migrate_text(text).splitlines()[:len(header)] == header


def test_example_valid():
    """The new example config passes the schema."""
    from snakemake.utils import validate
    validate(yaml.safe_load(load(EXAMPLE)), os.path.join(REPO, "schema", "config.schema.yaml"))


@pytest.mark.parametrize("scenario", OLD_SCENARIOS, ids=os.path.basename)
def test_scenario_merged(scenario):
    """Migrating the example and a scenario apart and merging them (as the
    dry runs do) is the same as migrating them merged."""
    ex, sc = load(OLD_EXAMPLE), load(scenario)
    together = merge(yaml.safe_load(ex), yaml.safe_load(sc) or {})
    check_round_trip(yaml.safe_dump(together, sort_keys=False))
    assert merge(migrated(ex), migrated(sc)) == migrated(
        yaml.safe_dump(together, sort_keys=False))


@pytest.mark.parametrize("step", [
    'STARsolo_kb_solo', 'STARsolo_cellsnp', 'STARsolo_demux', 'STARsolo_gt_demux',
    'STARsolo_split_bams_gt_demux', 'STARsolo_identify_swaps', 'STARsolo_gt_demux_multi_vcf',
    'all_multi_vcf', 'cellranger_count', 'cellranger_cellsnp', 'cellranger_demux',
    'cellranger_split_bams', 'multiome_alignment', 'multiome_gt_demux',
    'multiome_split_bams_gt_demux',
])
def test_last_step(step):
    text = f"last_step: {step}\n"
    check_derived(yaml.safe_load(text), migrated(text))


def test_comments_and_anchors():
    out = migrate_text(load(OLD_EXAMPLE))
    assert "whitelist: *3v3" in out
    assert "&gencode38o150" in out          # anchors no alias uses are kept
    assert "    # sjdboverhang: 89\n" in out  # comments follow their keys' indent
    assert "picard_metrics: all  # can be null, all, rnaseq or gc" in out


@pytest.mark.parametrize("text, error", [
    ("workflow:\n  aligner: STARsolo\n", "Already in the new layout"),
    ("last_step: STARsolo_demux\nnot_a_key: 1\n", "Unrecognised keys: not_a_key"),
    ("last_step: kallisto\n", "Unrecognised last_step"),
    ("demux_type: hashsolo\n", "Unrecognised demux_type"),
    ("a:\n  k: 1\n  k: 2\n", "Duplicate key"),
])
def test_errors(text, error):
    with pytest.raises(ValueError, match=error):
        migrate_text(text)
