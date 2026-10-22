#!/usr/bin/env python3

from typing import Any, Union # Need verion > 3.10
import pandas as pd, os


def get_demux_type(config) -> Union[None, str]:
    """Demultiplexing of this run: solo, vireo or both, add_solo or add_vireo
    (the method added to a run of the other tool, 'add_to_previous'), or
    None when no demux is needed."""
    d = config['experiment']['demultiplexing']
    method = d.get('method')
    if not method or str(method).lower() in ('none', 'null'):
        return None
    method = str(method).lower()
    return f"add_{method}" if d.get('add_to_previous') else method


def is_genotyped(config) -> bool:
    """Whether donor genotypes are used (else vireo names donors donor0, donor1, ...)."""
    return bool(config['experiment']['demultiplexing']['genotype']['enabled'])


# Output dirs of a demultiplexing run; 'run_tag' adds a sub-folder to each
RUN_TAGGED_DIRS = {
    'hashsolo': [
        'calico_solo_dir', 'final_count_matrix_dir', 'demultiplex_info_dir',
    ],
    'gt_demux': [
        'inp_for_cellsnp_dir', 'cellsnp_dir', 'vireosnp_dir',
        'final_count_matrix_dir', 'demultiplex_info_dir',
    ],
}


def apply_run_tag(config) -> None:
    """Add the 'run_tag' sub-folder to the dirs in RUN_TAGGED_DIRS.

    The untagged dirs are kept in config['_untagged_dirs'], so the outputs
    of other runs can still be found (see run_dir). Calling it again
    doesn't add the tag twice.
    """
    if '_untagged_dirs' in config:
        return

    config['_untagged_dirs'] = {
        section: {key: config['modules'][section][key] for key in keys}
        for section, keys in RUN_TAGGED_DIRS.items()
    }
    for section, keys in RUN_TAGGED_DIRS.items():
        for key in keys:
            config['modules'][section][key] = run_dir(
                config, section, key, config['experiment']['demultiplexing'].get('run_tag')
            )


# Dirs of the split bams pipeline; bams split for genotype checks
# (gt_check) and the final per-donor bams each get a sub-folder
SPLIT_PURPOSE_DIRS = [
    'inp_split_bams_dir', 'split_bams_dir', 'per_donor_split_log_dir',
    'new_temp_dir',
]


def split_purpose(config) -> str:
    return 'gt_check' if config['workflow']['gt_check'] else 'final'


def apply_split_purpose(config) -> None:
    """Add the split_purpose sub-folder to the dirs in SPLIT_PURPOSE_DIRS."""
    if '_split_purpose' in config:
        return
    config['_split_purpose'] = split_purpose(config)
    sb = config['modules']['split_bams']
    for key in SPLIT_PURPOSE_DIRS:
        sb[key] = f"{sb[key]}{config['_split_purpose']}/"


def run_dir(config, section: str, key: str, run_tag: Union[None, str]) -> str:
    """config['modules'][section][key] of the run tagged 'run_tag' (None: untagged)."""
    base = config['_untagged_dirs'][section][key]
    return f"{base}{run_tag}/" if run_tag else base


# def read_files_ext(fname, h: Union[None, list]=None) -> pd.DataFrame :
def read_files_ext(fname: str, **kwargs: Any) -> pd.DataFrame :    
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 

    if fname.endswith('.csv'):
        return pd.read_csv(fname, **kwargs)
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t', **kwargs)
    elif fname.endswith('.txt'):
        return pd.read_csv(fname, sep=' ', **kwargs)
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")


def ret_cols(fname: str) -> int :    
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 
    if fname.endswith('.csv'):
        return pd.read_csv(fname, nrows=1, header=None).shape[1]
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t', nrows=1, header=None).shape[1]
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")
    
