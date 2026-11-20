#!/usr/bin/env python3
"""Run one of WASP's python scripts with the numpy names it uses that newer
numpy removed restored.

WASP (github.com/bmvdgeijn/WASP, mapping/) was written for numpy 1.x. Of the
numpy names its scripts use, two are gone (else "module 'numpy' has no
attribute ..."):

    np.int   snptable.py; removed in numpy 1.24: the builtin int
    np.in1d  find_intersecting_snps.py; removed in numpy 2.4: np.isin of the
             flattened array, as np.in1d was

The rest of what they use (np.array, np.unique, ...; pysam's old names
Samfile, qname, pos, qual; PyTables' open_file, get_node) is still there,
and Python 3.12 runs them as they are: run end to end on the same inputs
with numpy 2.5 and 1.26 (python 3.12, pysam 0.24), the outputs are the same.
Only names missing from the numpy in use are set. Used by
helper_sh_scripts/wasp_donor.sh; any arguments after the script are its own:

    python helper_py_scripts/wasp_run.py $WASP_HOME/mapping/find_intersecting_snps.py \\
        --is_sorted --snp_index ... donor.bam
"""
import os
import runpy
import sys
import warnings

import numpy as np


def in1d(ar1, ar2, assume_unique=False, invert=False, **kwargs):
    """np.in1d as numpy < 2.4 had it: a 1-D result for the flattened ar1."""
    return np.isin(np.asarray(ar1).ravel(), ar2, assume_unique=assume_unique,
                   invert=invert, **kwargs)


# The removed names WASP's scripts use -> what they were
REMOVED = {'int': int, 'in1d': in1d}


def restore_aliases():
    """Set each removed name the numpy in use lacks; returns the names set.
    Where numpy still has them (np.int before 1.24, np.in1d before 2.4), its
    deprecation warnings are silenced, to keep them out of WASP's logs."""
    restored = []
    for name, value in REMOVED.items():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            present = hasattr(np, name)
        if not present:
            setattr(np, name, value)
            restored.append(name)
    warnings.filterwarnings('ignore', category=DeprecationWarning,
                            message=r'.*(np\.int\b|in1d).*')
    warnings.filterwarnings('ignore', category=FutureWarning,
                            message=r'.*(np\.int\b|in1d).*')
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
