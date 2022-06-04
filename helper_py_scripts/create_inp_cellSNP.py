#!/usr/bin/env python

import anndata as ad, string
import scanpy as sc, pandas as pd, numpy as np
import glob2, os, re, argparse
from collections import Counter
from openpyxl import load_workbook
from collections import defaultdict
# using datetime module
import datetime
import functools
from itertools import repeat


#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Produce inputs (text files) per sample for cellSNP, containing filtered barcodes for each")

parser.add_argument('inp', help="Path to the cached output of the final matrix (h5ad)")

# Optional parameters
parser.add_argument('-c', '--column', help="Column Name containing classifications of Doublets and Negatives. Default: SubID_cs", default='SubID_cs')
parser.add_argument('-o', '--output', help="Name of the output file")
parser.add_argument('--keep_all_cells', action='store_true', help="Use this flag when you want to retain cells not classified into a donor (e.g. doublets, negatives, etc.)")
parser.add_argument('-e', '--extra', nargs='+', help="Classficiation used when a cell is not present in the classification by hashsolo. NOTE: Accepts multi-word - space separated input. \
    If not used, use 'None'. Default: Not Present.", default=["Not", "Present"])
parser.add_argument('-d', '--doublet', help="Classficiation used when a cell is a doublet in the classification by hashsolo. If not used, use 'None'. Default: Doublet", default='Doublet')
parser.add_argument('-n', '--negative', help="Classficiation used when a cell is a negative in the classification by hashsolo. If not used, use 'None'. Default: Negative", default='Negative')
parser.add_argument('-b', '--barcode_len', type=int, help="Barcode length. For 10x it is 16", default=16)


args = parser.parse_args()

# Join
extra_val=' '.join(args.extra)


# Inform the user when the flag 'keep_all_cells' is used/turned on
if args.keep_all_cells:
    print("Flag 'keep_all_cells' is used. Hence, all cells will be considered!")


try:
    adata=ad.read(args.inp)
except:
    raise ValueError('The cached output(h5ad) of hashsolo doesn\'t exist!')
        
# Check if the given column exists in the h5ad file
if args.column.lower() not in map(str.lower, adata.obs_keys()):
    raise ValueError("Given value for the argument 'column' doesn't exist in the given h5ad file!")

# No. of Cells
cells=adata.n_obs

# Santiy check for the arguments
if not args.keep_all_cells:
    if extra_val not None and extra_val != "None":
        if extra_val.lower() not in adata.obs[args.column].str.lower():
            raise ValueError("The word(s) provided to the argument 'extra' is not found in the given column (argument 'column')")
        else:
            extra_val_ser=adata.obs[args.column].str.lower() != extra_val.lower()
    else:
        extra_val_ser=pd.Series(repeat(True, cells))

    if args.doublet.lower() not None and args.doublet.lower() != "None":
        if args.doublet.lower() not in adata.obs[args.column].str.lower():
            raise ValueError("The word(s) provided to the argument 'doublet' is not found in the given column (argument 'column')")
        else:
            doublet_ser=adata.obs[args.column].str.lower() != args.doublet.lower()
    else:
        doublet_ser=pd.Series(repeat(True, cells))

    if args.negative.lower() not None and args.negative.lower() != "None":
        if args.negative.lower() not in adata.obs[args.column].str.lower():
            raise ValueError("The word(s) provided to the argument 'negative' is not found in the given column (argument 'column')")
        else:
            negative_ser=adata.obs[args.column].str.lower() != args.negative.lower()
    else:
        negative_ser=pd.Series(repeat(True, cells))


    adata=adata[ extra_val_ser & doublet_ser & negative_ser ].copy()



with open(args.output, 'w') as f:
    for item in adata.obs_names.tolist():
        item=item[:args.barcode_len] # To remove the suffixes used to make barcodes unique, if present
        f.write("%s\n" % item)