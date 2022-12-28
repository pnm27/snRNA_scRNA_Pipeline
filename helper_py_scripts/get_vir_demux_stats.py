#!/usr/bin/env python3

import pandas as pd
import os, yaml
import glob2
from time import sleep
from typing import Union


"""This script is to create a file with vir outputs for a list of samples.

This script serves as an intermediate or temporary use to delineate the demultiplexing
results of vireo run on a list of samples, which are present inside a text file or a yaml file.
Usage:
    python3 ./get_vir_demux_stats.py

Author:
    Prashant N M
"""

def read_inp(x: str):
    """ This function identifies the type of the input file and returns a pandas dataframe
    """

    if x.endswith('.yaml') or x.endswith('.yml'):
        with open(x) as fout:
            sample_dict = yaml.load(fout, Loader=yaml.SafeLoader)

        if isinstance(sample_dict, list):
            print(f"Samples contained in the input file: {len(sample_dict)}")
            return pd.DataFrame(sample_dict, columns=["col1"])

        elif isinstance(sample_dict, dict):
            temp_val=[]
            t=0
            for k,v in sample_dict.items():
                if isinstance(v, list):
                    t+=1
                    if t>1:
                        print("Multiple sets of samples present. They will be consolidated into one")

                    temp_val.extend(v)
                elif isinstance(v, str):
                    temp_val.append(v)
                else:
                    raise ValueError("Input yaml file contains at least 2 levels of dict! Max allowed is 1")

            print(f"Samples contained in the input file: {len(temp_val)}")
            return pd.DataFrame(temp_val, columns=["col1"])

    elif x.endswith('.txt') or os.path.isfile(x):
        f = pd.read_csv(x, names=["col1"])
        print(f"Samples contained in the input file: {f.shape[0]}")
        return f

    else:
        raise ValueError("Wrong input file provided! Provide either a yaml file or a text file (no header)")


# Project-specific
def conv_type(x):
    try:
        return x.strip()
    except:
        return x


def ret_strings(x, sep=','):
    return sep.join(x)


def get_don_ids(x, t_df) -> Union[str, int]:
    if x == 'Doublet' or x == 'Negative' or x.startswith('donor') or t_df is None:
        return x
    else:
        return t_df.loc[t_df.iloc[:, 1:].isin([str(x).strip()]).any(axis=1), "SubID"].values[0]


# The input file used to run Snakemake
all_files = read_inp("/sc/arion/projects/CommonMind/pnm/choroid_plexus/fastq_files.txt")
vir_name_conv = False # Whether to convert the donor names in vir output or not
conv_file = None # Converter file, if donor names in the vir output needs to be converted
vir_out_dir = "/sc/arion/projects/CommonMind/pnm/choroid_plexus/demux_gt/vireoSNP_1kGP_isec_rerun/"
out_file = "/sc/arion/projects/CommonMind/pnm/choroid_plexus/Vir_demux_stats_round1.txt"

if vir_name_conv and conv_file is not None:
    # Just to be safe convert everything to string
    conv_df = pd.read_csv(conv_file, dtype=str)
    all_cols=list(conv_df)
    conv_df[all_cols] = conv_df[all_cols].applymap(conv_type)

a = []
for v in all_files["col1"]:
    # Sample name processing, if required
    vals=v.split("/")
    fname = vals[-1].replace('-cDNA', '')
    
    temp_f = pd.read_csv(os.path.join(vir_out_dir, fname) + "/summary.tsv", sep='\t')
        
    t_dons_l = temp_f["Var1"].tolist()
    t_cc_l = temp_f["Freq"].tolist()
    dons=[]
    cells=[]
    doub_val=0
    neg_val=0
    for i, j in enumerate(t_dons_l):
        if j != 'doublet' and j != 'unassigned':
            dons.append(j)
            cells.append(str(t_cc_l[i]))
        elif j == 'doublet':
            doub_val = t_cc_l[i]
        elif j == 'unassigned':
            neg_val = t_cc_l[i]

    if vir_name_conv and conv_file is not None:
        sid=[]
        for d in dons:
            if not d.startswith('donor'):
                sid.append(get_don_ids('_'.join(d.split('_')[1:]), conv_df))
            else:
                sid.append(d)

if vir_name_conv and conv_file is not None:
    a.append((fname, ret_strings(dons), ret_strings(sid), ret_strings(cells), neg_val, doub_val ))
    a_df = pd.DataFrame(a, columns=["sample", "donors", "SubIDs", "cell_counts", "Negatives", "Doublets"])
else:
    a.append((fname, ret_strings(dons), ret_strings(cells), neg_val, doub_val ))
    a_df = pd.DataFrame(a, columns=["sample", "donors", "cell_counts", "Negatives", "Doublets"])

a_df.to_csv(out_file, sep = " ", index=False)