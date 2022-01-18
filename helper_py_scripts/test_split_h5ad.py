#!/usr/bin/env python

import re, glob2, os, argparse, gc
import pandas as pd, anndata as ad, scanpy as sc
from pathlib import Path
from collections import Counter

samples = sorted(glob2.glob("/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/h5ad/*.h5ad"))

for sample in samples:
	try:
		adata = ad.read(sample)
        
	except:
		print("Problem loading donor: {}".format(sample))
		os.remove(sample)
