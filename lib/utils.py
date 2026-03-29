#!/usr/bin/env python3

from typing import Any, Union # Need verion > 3.10
import pandas as pd, os


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
    
