#!/usr/bin/env python3
"""Run one of WASP's python scripts with numpy's removed type aliases restored.

WASP (github.com/bmvdgeijn/WASP) uses np.int (mapping/snptable.py), which
numpy 1.24 removed, so with a newer numpy find_intersecting_snps.py fails
with "module 'numpy' has no attribute 'int'". The aliases were the builtins
(np.int is int), so restoring them changes nothing else. Used by
helper_sh_scripts/wasp_donor.sh; any arguments after the script are its own:

    python helper_py_scripts/wasp_run.py $WASP_HOME/mapping/find_intersecting_snps.py \\
        --is_sorted --snp_index ... donor.bam
"""
import os
import runpy
import sys

import numpy as np

# numpy < 1.24 has them (deprecated); with numpy >= 1.24 they raise AttributeError
ALIASES = {'int': int, 'float': float, 'bool': bool, 'object': object,
           'complex': complex, 'str': str}


def restore_aliases():
    """Set each removed alias; returns the names set."""
    restored = []
    for name, builtin in ALIASES.items():
        try:
            getattr(np, name)
        except AttributeError:
            setattr(np, name, builtin)
            restored.append(name)
    return restored


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: wasp_run.py <WASP script> [its arguments]")
    script = sys.argv[1]
    if not os.path.isfile(script):
        sys.exit(f"wasp_run.py: no script {script}")
    restore_aliases()
    # As 'python script args': its own argv, and its folder first on the path
    # (find_intersecting_snps.py imports snptable, util from there)
    sys.argv = sys.argv[1:]
    sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
    runpy.run_path(script, run_name='__main__')


if __name__ == '__main__':
    main()
