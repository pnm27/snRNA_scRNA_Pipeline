#!/usr/bin/env python

import scanpy as sc, anndata as ad, pandas as pd, numpy as np
import glob2, os, re

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()

a = glob2.glob("/sc/arion/projects/psychAD/HTO_info/feature_barcodes/round*/*/bus_count_out/output.mtx")
for i in range(len(a)):
   genes=sc.read_mtx(a[i])
   genes.var_names=pd.read_csv(a[i].replace('.mtx', '.genes.txt'), header=None)[0]
   genes.obs_names=pd.read_csv(a[i].replace('.mtx', '.barcodes.txt'), header=None)[0]
   r = re.search('/round([0-9]+)/Sample', a[i]).group(1)
   s = re.search('/Sample_(NPSAD-.*)/bus', a[i]).group(1)
   genes.write('cache/round' + r + '_Sample-' + s + '.h5ad')

