"""Defaults for a config of existing bams ('workflow: start: split_bams' or
'wasp').

Such a run reads only the split_bams (or wasp) settings, but the schema
needs every section. So its config goes over new_config.example.yaml and
the defaults of its start (schema/split_bams_defaults.yaml,
schema/wasp_defaults.yaml), and needs only what differs: e.g. the external
labels and, for bams outside the aligner's folders, the pools table; or
the WASP samples table, genotypes and remap genome.
Configs of runs from the fastqs ('start: align') are left as they are, so a
missing genome or fastq dir is an error, not a default.
"""

import copy
import os
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(REPO, "new_config.example.yaml")
SPLIT_BAMS_DEFAULTS = os.path.join(REPO, "schema", "split_bams_defaults.yaml")
WASP_DEFAULTS = os.path.join(REPO, "schema", "wasp_defaults.yaml")
# 'workflow: start' -> its defaults
START_DEFAULTS = {'split_bams': SPLIT_BAMS_DEFAULTS, 'wasp': WASP_DEFAULTS}


def merge(base: dict, over: dict) -> dict:
    """'over' on top of 'base': nested dicts merged, other values replaced."""
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _start(config: dict):
    workflow = config.get('workflow')
    return workflow.get('start') if isinstance(workflow, dict) else None


def starts_at_split_bams(config: dict) -> bool:
    return _start(config) == 'split_bams'


def start_defaults(start: str) -> dict:
    """The example config with the defaults of a 'start' on top."""
    with open(EXAMPLE) as f:
        example = yaml.safe_load(f)
    with open(START_DEFAULTS[start]) as f:
        return merge(example, yaml.safe_load(f))


def split_bams_defaults() -> dict:
    return start_defaults('split_bams')


def fill_defaults(config: dict) -> dict:
    """Fill in (in place) the keys a 'start: split_bams' (or 'wasp') config
    leaves out."""
    start = _start(config)
    if start in START_DEFAULTS:
        filled = merge(start_defaults(start), config)
        config.clear()
        config.update(filled)
    return config
