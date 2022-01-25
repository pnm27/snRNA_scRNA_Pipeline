#!/usr/bin/env python

import scanpy as sc, anndata as ad, pandas as pd, numpy as np
import glob2, os, re, argparse

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()

parser = argparse.ArgumentParser(description="Split a(set of) h5ad(s) into multiple h5ad based on a obs column in it(them)")

parser.add_argument('matrix_file', help="Path to output.mtx")
parser.add_argument('genes_file', help="Path to output.genes.txt")
parser.add_argument('barcodes_file', help="Path to output output.barcodes.txt")

# Optional parameters
parser.add_argument('-o', '--output', help="Name of the output file (h5ad)")

args = parser.parse_args()

inp_mat = args.matrix_file
inp_genes = args.genes_file
inp_bc = args.barcodes_file
if args.output.endswith('.h5ad'):
	out_file=args.output
else:
	out_file=args.output + '.h5ad'

genes=sc.read_mtx(inp_mat)
genes.var_names=pd.read_csv(inp_genes, header=None)[0]
genes.obs_names=pd.read_csv(inp_bc, header=None)[0]
genes.write(out_file)

