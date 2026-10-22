#!/usr/bin/env python3
"""Check the config against schema/config.schema.yaml and explain what's wrong.

snakemake's validate() stops at the first error, with jsonschema's message.
check_config() lists every error, each with the reason from the schema: the
'x-explain' of the rule it breaks, or the key's 'description'. The rules
across keys are in the schema's 'allOf'; its 'x-warnings' are only printed.
It runs before validate(), which then fills in the defaults.

Also, duplicate keys in a config file: YAML keeps the last value, so the
first one is silently ignored.
"""

import sys
import yaml
from jsonschema import validators


class ConfigError(ValueError):
    pass


def _where(path) -> str:
    return ': '.join(str(p) for p in path) or '(top level)'


def _reason(schema, schema_path) -> str:
    """The 'x-explain' (or 'description') nearest to the failing keyword."""
    node, reason = schema, None
    for key in schema_path:
        try:
            node = node[key]
        except (KeyError, IndexError, TypeError):
            break
        if isinstance(node, dict):
            reason = node.get('x-explain') or node.get('description') or reason
    return reason


def _message(error, schema) -> str:
    if error.validator == 'additionalProperties':
        known = sorted(error.schema.get('properties', {}))
        unknown = sorted(set(error.instance) - set(known))
        return (
            f"Unknown key(s) {', '.join(repr(k) for k in unknown)} in "
            f"{_where(error.path)}.\n  Why: only these keys are read there: "
            f"{', '.join(known)}. A typo, or a key of the old layout? "
            "(helper_py_scripts/migrate_config.py converts old configs)"
        )
    if error.validator == 'not':
        msg = f"{_where(error.path)}: {error.instance!r} isn't allowed here"
    else:
        msg = f"{_where(error.path)}: {error.message}"
    reason = _reason(schema, error.schema_path)
    return f"{msg}\n  Why: {reason}" if reason else msg


def config_problems(config, schema):
    """(errors, warnings) of 'config', as explained messages."""
    cls = validators.validator_for(schema)
    errors = [
        _message(e, schema)
        for e in sorted(cls(schema).iter_errors(config), key=lambda e: list(e.path))
    ]
    warnings = []
    for rule in schema.get('x-warnings', []):
        cond = {k: v for k, v in rule.items() if k in ('if', 'then')}
        if not cls(cond).is_valid(config):
            warnings.append(rule['x-warn'])
    return errors, warnings


class _UniqueKeyLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader, node, deep=False):
    seen = {}
    for key_node, _ in node.value:
        if key_node.tag == 'tag:yaml.org,2002:merge':  # '<<: *anchor'
            continue
        key = loader.construct_object(key_node, deep=deep)
        if key in seen:
            raise ConfigError(
                f"Duplicate key '{key}' (lines {seen[key] + 1} and "
                f"{key_node.start_mark.line + 1}).\n  Why: YAML keeps the "
                "last value, so the first one is silently ignored; remove one"
            )
        seen[key] = key_node.start_mark.line
    return loader.construct_mapping(node, deep=deep)


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def duplicate_keys(path):
    """Message for the first duplicate key in a YAML file, if any."""
    try:
        with open(path) as f:
            yaml.load(f, Loader=_UniqueKeyLoader)
    except ConfigError as e:
        return f"{path}: {e}"
    return None


def check_config(config, schema_path, config_files=()):
    """Raise ConfigError listing every problem; print the warnings."""
    with open(schema_path) as f:
        schema = yaml.safe_load(f)
    errors = [m for m in map(duplicate_keys, config_files) if m]
    if 'last_step' in config or 'workflow' not in config:
        errors.append(
            "The config is in the old layout.\n  Why: the pipeline reads "
            "the new layout (references, workflow, experiment, modules, "
            "execution). Convert it with: python3 "
            "helper_py_scripts/migrate_config.py <config> -o <new config>"
        )
    else:
        problems, warnings = config_problems(config, schema)
        errors += problems
        for w in warnings:
            print(f"Config warning: {w}", file=sys.stderr)
    if errors:
        raise ConfigError(
            f"{len(errors)} problem(s) in the config "
            "(see docs/source/config_caveats.md):\n\n" + "\n\n".join(errors)
        )
