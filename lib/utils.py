#!/usr/bin/env python3

from typing import Any, Union # Need verion > 3.10
import pandas as pd, os


def get_demux_type(config) -> Union[None, str]:
    """Return config['demux_type'] lowercased, or None when no demux is needed."""
    demux_type = config.get('demux_type')
    if not demux_type or str(demux_type).lower() in ('none', 'null'):
        return None
    return str(demux_type).lower()


# Output dirs of a demultiplexing run; 'run_tag' adds a sub-folder to each
RUN_TAGGED_DIRS = {
    'hashsolo_demux_pipeline': [
        'calico_solo_dir', 'final_count_matrix_dir', 'demultiplex_info_dir',
    ],
    'gt_demux_pipeline': [
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
        section: {key: config[section][key] for key in keys}
        for section, keys in RUN_TAGGED_DIRS.items()
    }
    for section, keys in RUN_TAGGED_DIRS.items():
        for key in keys:
            config[section][key] = run_dir(
                config, section, key, config.get('run_tag')
            )


def run_dir(config, section: str, key: str, run_tag: Union[None, str]) -> str:
    """config[section][key] of the run tagged 'run_tag' (None: untagged)."""
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
    
