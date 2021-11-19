#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

# Solo didn't run through scvi, scvi-tools nor scanpy.external
# Only this seems to work
from solo import hashsolo
import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import glob2, os, re, argparse
from collections import Counter
from openpyxl import load_workbook
from collections import defaultdict


sc.settings.set_figure_params(dpi_save=400, format='png', color_map = 'viridis_r')
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()

parser = argparse.ArgumentParser(description="Create h5ad output after running calico_solo(hashsolo)")

parser.add_argument('bustools_out', help="Path to cached output of bustools kite pipeline(h5ad)")
parser.add_argument('matrix_file', help="Path to matrix.mtx.gz")
parser.add_argument('output_file', help="Path to store output file(h5ad)")

# Optional parameters
parser.add_argument('-m', '--max_mito', type=int, help="Max mitochondrial genes(in percent) per cell. Default: 5", default=5)
parser.add_argument('-g', '--min_genes', type=int, help="Min #genes per cell. Default: 1000", default=1000)
parser.add_argument('-c', '--min_cells', type=int, help="Min #cells expressing a gene for it to pass the filter. Default: 10", default=10)


args = parser.parse_args()


# Parameters for filtering
max_mito = args.max_mito
min_genes = args.min_genes
min_cells = args.min_cells


t2g = pd.read_csv("/sc/arion/projects/psychAD/pnm/Hs_allchr_MT.txt", skiprows=1, usecols=range(2),names=["gene_id", "gene_name"], sep="\t")
t2g.index = t2g.gene_id
t2g = t2g.loc[~t2g.index.duplicated(keep='first')]

hash_data = ad.read(args.bustools_out)
adata = sc.read_10x_mtx(args.matrix_file[:-13], make_unique=True, var_names= "gene_ids", cache=True)
adata.var_names_make_unique()
adata.var["gene_id"] = adata.var.index.values
adata.var["gene_name"] = adata.var.gene_id.map(t2g["gene_name"])
adata.var_names = adata.var_names.to_series().map(lambda x: x + '_index')
adata = adata[:, pd.notna(adata.var["gene_name"])]
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')
adata.obs_names = adata.obs_names.to_series().map(lambda x: re.sub('-.*', '', x))
   
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
# Filter data wrt mito content
adata.var["mito"] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
adata = adata[adata.obs["pct_counts_mito"]< max_mito, :]

# Some barcodes in gene expression matrix doesn't exist in the bus output
try:
    hash_data = hash_data[adata.obs_names].copy()
except:
    filt_bc = [b for b in adata.obs_names.to_list() if b in hash_data.obs_names]
    hash_data = hash_data[filt_bc].copy()

hashsolo.hashsolo(hash_data)
hash_data.write(args.output_file)
'''

Default:
priors: list = [.01, .8, .19],
pre_existing_clusters: str = None,
clustering_data: anndata.AnnData = None,
resolutions: list = [.1, .25, .5, .75, 1],
number_of_noise_barcodes: int = None,
inplace: bool = True

'''
