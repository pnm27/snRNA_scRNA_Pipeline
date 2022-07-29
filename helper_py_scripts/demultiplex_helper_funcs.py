#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

import pandas as pd, numpy as np, os



# Basic helper functions-------------------------------------------------------------------------------------------
def auto_read(fname) -> pd.DataFrame :
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 
    if fname.endswith('.csv'):
        return pd.read_csv(fname)
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t')
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")


# Simplify these 'parse' functions
def parse_HTO(wet_lab_df, col_val, fname, s_name, hs=None) -> list:
    if isinstance(col_val, str):
        htos = wet_lab_df[col_val]
        test_len = len(htos)
        # No command-line params and inference that all HTOs are present in one row separated by ","
        if hs == None and test_len == 1 and htos.str.count(',').values[0] > 1:
            return htos.values[0].split(',')
        # No command-line params and inference that all HTOs are present in one row separated by whitespaces    
        elif hs == None and test_len == 1 and len(htos.split()) > 1:
            return htos.values[0].split(',')
        elif hs == None and test_len > 1:
            return htos.tolist()
        elif hs != None and test_len == 1 and htos.str.count(hs).values[0] > 1:
            return htos.values[0].split(hs)
        elif hs != None and test_len > 1:
            raise ValueError(f"After subsetting sample {s_name} from the wet lab file {fname}, there are multiple rows ({test_len}) of HTOs for this sample while a separator value is also provided.")
        elif hs != None and test_len == 1 and htos.str.count(hs).values[0] == 1:
            raise ValueError(f"Either the given separator {hs} is wrong or the sample {s_name} has incomplete HTO values in the wet lab file {fname}")
        else:
            raise ValueError(f"Something is wrong with the given input(s):\n\twet lab file: {fname}\n\tsample: {s_name}\n\tHTO-separator: {hs}")

    # If 'multiple' HTOs are used for identification
    elif isinstance(col_val, list):
        # If col_val is a list of headers
        try:
            return wet_lab_df[col_val].to_numpy().tolist()
        # If col_val is a list of integers representing columns
        except:
            return wet_lab_df.iloc[:, [col_val]].to_numpy().tolist()


    else:
        raise ValueError(f"Unexpected typing for the columns! Expected either str or list but got {type(col_val)}")


# Assumes similar construct like the HTOs
def parse_subids(wet_lab_df, col_val, fname, s_name, hs=None) -> list:
    sub = wet_lab_df[col_val]
    test_len = len(sub)
    # No command-line params and inference that all HTOs are present in one row separated by ","
    if hs == None and test_len == 1 and sub.str.count(',').values[0] > 1:
        return sub.values[0].split(',')
    # No command-line params and inference that all HTOs are present in one row separated by whitespaces    
    elif hs == None and test_len == 1 and len(sub.split()) > 1:
        return sub.values[0].split(',')
    elif hs == None and test_len > 1:
        return sub.tolist()
    elif hs != None and test_len == 1 and sub.str.count(hs).values[0] > 1:
        return sub.values[0].split(hs)
    elif hs != None and test_len > 1:
        raise ValueError(f"After subsetting sample {s_name} from the wet lab file {fname}, there are multiple rows ({test_len}) of HTOs for this sample while a separator value is also provided.")
    elif hs != None and test_len == 1 and sub.str.count(hs).values[0] == 1:
        raise ValueError(f"Either the given separator {hs} is wrong or the sample {s_name} has incomplete donor ids in the wet lab file {fname}")
    else:
        raise ValueError(f"Something is wrong with the given input(s):\n\twet lab file: {fname}\n\tsample: {s_name}\n\tHTO-separator: {hs}")