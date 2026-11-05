"""Defaults for a config of existing bams ('workflow: start: split_bams').

Such a run reads only the split_bams settings, but the schema needs every
section. So its config goes over new_config.example.yaml and
schema/split_bams_defaults.yaml, and needs only what differs: the external
labels and, for bams outside the aligner's folders, the pools table.
Configs of runs from the fastqs ('start: align') are left as they are, so a
missing genome or fastq dir is an error, not a default.
"""

import copy
import os
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(REPO, "new_config.example.yaml")
SPLIT_BAMS_DEFAULTS = os.path.join(REPO, "schema", "split_bams_defaults.yaml")


def merge(base: dict, over: dict) -> dict:
    """'over' on top of 'base': nested dicts merged, other values replaced."""
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def starts_at_split_bams(config: dict) -> bool:
    workflow = config.get('workflow')
    return isinstance(workflow, dict) and workflow.get('start') == 'split_bams'


def split_bams_defaults() -> dict:
    """The example config with the defaults of a split_bams run on top."""
    with open(EXAMPLE) as f:
        example = yaml.safe_load(f)
    with open(SPLIT_BAMS_DEFAULTS) as f:
        return merge(example, yaml.safe_load(f))


def fill_defaults(config: dict) -> dict:
    """Fill in (in place) the keys a 'start: split_bams' config leaves out."""
    if starts_at_split_bams(config):
        filled = merge(split_bams_defaults(), config)
        config.clear()
        config.update(filled)
    return config
