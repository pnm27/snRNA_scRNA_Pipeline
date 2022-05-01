#!/usr/bin/env python

import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import glob2, os, re, argparse
from collections import Counter
from openpyxl import load_workbook
from collections import defaultdict
# using datetime module
import datetime
import functools


#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Produce inputs (text files) per sample for cellSNP, containing filtered barcodes for each")

parser.add_argument('inp', help="Path to the cached output of the final matrix (h5ad)")

# Optional parameters
parser.add_argument('-c', '--column', help="Column Name containing classifications of Doublets and Negatives. Default: SubID_cs", default='SubID_cs')
parser.add_argument('-o', '--output', help="Name of the output file")
parser.add_argument('-e', '--extra', nargs='+', help="Classficiation used when a cell is not present in the classification by hashsolo. NOTE: Accepts multi-word, space separated input. Default: Not Present", default=["Not", "Present"])
parser.add_argument('-d', '--doublet', help="Classficiation used when a cell is a doublet in the classification by hashsolo. Default: Doublet", default='Doublet')
parser.add_argument('-n', '--negative', help="Classficiation used when a cell is a negative in the classification by hashsolo. Default: Negative", default='Negative')
parser.add_argument('-b', '--barcode_len', type=int, help="Barcode length. For 10x it is 16", default=16)

args = parser.parse_args()

extra_val=' '.join(args.extra)

try:
    adata=ad.read(args.inp)
except:
    raise ValueError('The cached output(h5ad) of hashsolo doesn\'t exist!')
        

adata=adata[(adata.obs[args.column].str.lower() != args.doublet.lower()) & (adata.obs[args.column].str.lower() != args.negative.lower()) & (adata.obs[args.column].str.lower() != extra_val.lower())].copy()
with open(args.output, 'w') as f:
    for item in adata.obs_names.tolist():
        item=item[:args.barcode_len] # To remove the suffixes used to make barcodes unique, if present
        f.write("%s\n" % item)