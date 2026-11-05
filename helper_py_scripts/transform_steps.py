#!/usr/bin/env python3
"""Transformation steps of the JSON/YAML recipes ('pythonic' and 'regex').

The steps are defined in schema/demul_samples_annotate_h5ad.schema.json
($defs: pythonic_step, regex_step), and copied into
schema/config.schema.yaml for 'modules: split_bams: split_by: external'.
Used by demultiplex_helper_funcs.process_columns (h5ad annotations) and
normalize_external_labels.py (barcodes and pool names).

Example: the part of 'P12_L1_AAACCTGAGC-1' after the last '_'::

    apply_steps('P12_L1_AAACCTGAGC-1', [
        {'type': 'pythonic', 'operation': 'split', 'by': '_'},
        {'type': 'pythonic', 'operation': 'select_index', 'index': -1},
    ])
"""

from __future__ import annotations
import re
from typing import Any, Iterable

try:  # 'regex' supports more than 're' (e.g. variable-length lookbehinds)
    import regex as _re
except ImportError:  # pragma: no cover
    _re = re

# Operations without arguments: str methods of the same name
_STR_METHODS = {'lower', 'upper', 'title', 'strip', 'rstrip', 'lstrip'}


def apply_pythonic(value: Any, t: dict[str, Any]) -> Any:
    r"""Apply one 'pythonic' step (see the module docstring).

    Parameters
    ----------
    value
        A string, or a list of them (after 'split')
    t
        The step: 'operation' and its arguments ('by', 'index', 'start',
        'end', 'old', 'new')

    Returns
    -------
    Any
        The transformed value
    """
    op = t["operation"]
    if op in _STR_METHODS:
        return getattr(value, op)()
    if op == "split":
        return value.split(t.get("by", " "))
    if op == "join":
        return t.get("by", " ").join(value)
    if op == "select_index":
        return value[t["index"]]
    if op == "select_range":
        return value[t.get("start"): t.get("end")]
    if op == "replace":
        return value.replace(t["old"], t["new"])
    raise ValueError(f"Unknown pythonic operation: '{op}'")


def _flags(flags: Any) -> int:
    """'IGNORECASE' or 'IGNORECASE|MULTILINE' (or an int) as regex flags."""
    if not flags:
        return 0
    if isinstance(flags, int):
        return flags
    out = 0
    for name in str(flags).split('|'):
        out |= getattr(_re, name.strip().upper())
    return out


def apply_regex(value: Any, t: dict[str, Any]) -> str:
    r"""Apply one 'regex' step: re.sub(pattern, replacement, value).

    Parameters
    ----------
    value
        The string to transform
    t
        The step: 'pattern', 'replacement' and, optionally, 'flags'
        (e.g. IGNORECASE)

    Returns
    -------
    str
        The transformed value
    """
    return _re.sub(t["pattern"], t["replacement"], value,
                   flags=_flags(t.get("flags")))


def apply_steps(value: Any, steps: Iterable[dict[str, Any]] | None) -> Any:
    """Apply 'steps' (pythonic or regex) in order."""
    for t in steps or []:
        if t["type"] == "pythonic":
            value = apply_pythonic(value, t)
        elif t["type"] == "regex":
            value = apply_regex(value, t)
        else:
            raise ValueError(f"Unknown transformation type: '{t['type']}'")
    return value
