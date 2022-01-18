#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

# Solo didn't run through scvi, scvi-tools nor scanpy.external
# Only this seems to work
from solo import hashsolo
import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import glob2, os, re
from collections import Counter
from openpyxl import load_workbook
from collections import defaultdict
import subprocess, shlex



sc.settings.set_figure_params(dpi_save=400, format='png', color_map = 'viridis_r')
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()

dir_path = snakemake.input[0].replace('/matrix.mtx.gz', '')
barcodes_file = os.path.join(dir_path, 'barcodes.tsv.gz')
samp_uniq = re.search('/Sample_(NPSAD-.*)/counts/', snakemake.input[0]).group(1)
comm_l = f"comm -12 <(cat {snakemake.params[0]}) <(zcat {barcodes_file}) > /sc/arion/projects/psychAD/pnm/{samp_uniq}_test_bc.tsv"
outp = subprocess.run(comm_l, shell=True, capture_output=True, executable="/bin/bash")

print(outp.returncode)

hash_data=ad.read_mtx(os.path.join(dir_path, 'matrix.mtx.gz'))
hash_data=hash_data.T
var_names=pd.read_csv(os.path.join(dir_path, 'features.tsv.gz'), compression='gzip', sep='\t', header=None)
obs_names=pd.read_csv(barcodes_file, compression='gzip', sep='\t', header=None)
red_bcs=pd.read_csv(f"/sc/arion/projects/psychAD/pnm/{samp_uniq}_test_bc.tsv",  sep='\t', header=None)
hash_data.var.index = var_names.iloc[:, 0].tolist()
hash_data.obs.index = obs_names.iloc[:, 0].tolist()
hash_data = hash_data[red_bcs.iloc[:, 0].tolist()].copy()

print(hash_data)
hash_data.write(snakemake.output[0])
#subprocess.run(f"rm /sc/arion/projects/psychAD/pnm/{samp_uniq}_test_bc.tsv", shell=True, capture_output=True, executable="/bin/bash")
