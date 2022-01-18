#!/sc/arion/work/prashf01/conda/envs/scanpy/bin/python


import scanpy as sc, anndata as ad, pandas as pd, numpy as np
import glob2

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()

#samples_file=pd.read_csv("/sc/arion/projects/psychAD/pnm/fastq_files.txt",names=["samples"], sep="\t")
#sample_paths = [glob2.glob('/sc/arion/projects/psychAD/STARsolo_bams/' + sample + '*_Solo.out/GeneFull/filtered_Lun/') for sample in samples_file.samples]
#samples = samples_file.samples.map(lambda x: x.split('/')[1].replace('-cDNA', ''))
sample_paths = glob2.glob('/sc/arion/projects/psychAD/Single_cell_data/alignment/*/cDNA/results/*/outs/filtered_feature_bc_matrix')
#[sc.read_10x_mtx(sample_path[0], make_unique=True, var_names= "gene_ids", cache=True) for sample_path in sample_paths]
[sc.read_10x_mtx(sample_path, make_unique=True, var_names= "gene_ids", cache=True) for sample_path in sample_paths]
