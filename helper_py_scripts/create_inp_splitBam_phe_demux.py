#!/usr/bin/env python
from itertools import repeat
import pandas as pd, numpy as np, anndata as ad, scanpy as sc
import sys, os, re, glob2, argparse
from collections import Counter


# Save dataframe to file
def save_df(df, suff, op):
	if suff.lower() == '.csv':
		df.to_csv(op, index=False)
	elif suff.lower() == '.tsv' or suff.lower() == '.txt':
		df.to_csv(op, sep = "\t", index=False)
	else:
		print("Can't save to file because of incorrect extension. Extension can be 'csv', 'tsv' or 'txt'")



vir_class = pd.read_csv(snakemake.input[0], sep='\t', usecols=["cell", "donor_id"])
vir_class = vir_class[(vir_class["donor_id"] != "doublet") & (vir_class["donor_id"] != "unassigned")]
vir_class.reset_index(drop=True, inplace=True)
# vir_class.rename(columns={"cell":"barcodes", "donor_id":"Subj_ID"}, inplace=True, errors="raise")
vir_class["donor_id"] = vir_class["donor_id"].apply(lambda x: '_'.join(x.split('_')[1:]) if not x.startswith('donor') else x)

# For converting genotype IDs to Subject IDs-------------------------------
conv_df = pd.read_csv(snakemake.params['conv_df'])
conv_df = conv_df.loc[conv_df['primary_genotype'].isin(vir_class['donor_id']), ["SubID", "primary_genotype"]]
vir_class['Subj_ID'] = vir_class['donor_id'].apply(lambda x: conv_df.loc[conv_df["primary_genotype"] == x, "SubID"].values[0] if not x.startswith('donor') else x)
del vir_class['donor_id']
vir_class.rename(columns={"cell":"barcodes"}, inplace=True, errors="raise")
# -------------------------------------------------------------------------
vir_class = vir_class[["Subj_ID", "barcodes"]]


file_ext = re.search(r'(\.[^.]+)$', snakemake.output[0]).group(1)


save_df(vir_class, file_ext, snakemake.output[0])