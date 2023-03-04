#!/usr/bin/env python3

"""Creates a list of effective barcodes for cellSNP run.

This script takes a count matrix as an input and creates a list of 
'effective' barcodes, dependent on the condition if some cells should 
be removed or not (e.g. if already a demultiplexing tool was run with an
h5ad column containing classifications and if some of the cells that
have been classified as a doublet and/or negative cells have to be removed 
then too this script can be used).

"""

import anndata as ad, string
import scanpy as sc, pandas as pd, numpy as np
import os, re, sys, argparse
from collections import Counter
from collections import defaultdict
# using datetime module
from itertools import repeat


def get_argument_parser():
    """Generate and return argument parser."""

    # Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Produce inputs "
    "(text files) per sample for cellSNP, containing filtered barcodes "
    "for each. If 'prev' flag is used then an h5ad file is expected"
    " and if not then a path containing the 10 mtx files. For an h5ad "
    "file one can specifiy if cell classifications are present or not "
    "and if present, then whether certain cell classes need to be "
    "removed or not ('keep_all_cells' flag)"
    )

    parser.add_argument('inp', help="Path to the cached output of the final "
    "matrix (h5ad) or Path containing 10x mtx files.",
    )
    parser.add_argument('-b', '--barcode_len', nargs='?', 
    type=int, help="Barcode length. For 10x, it is 16. When parameter "
    "present but no value provided: 16. Default: None", 
    const=16, default=None,
    )

    prev_count_file = parser.add_argument_group('h5ad_file')

    # Optional parameters
    prev_count_file.add_argument('--prev', action='store_true', 
    help="Use this flag when the input is an h5ad with cells "
    "classification present.",
    )
    prev_count_file.add_argument('-c', '--column', nargs='?', 
    help="Column name "
    "containing classifications of Doublets and Negatives. Value "
    "when parameter present but no value provided: SubID_cs. Default: None",
    const='SubID_cs', default=None,
    )
    prev_count_file.add_argument('-o', '--output', 
    help="Name of the output file",
    default="inp_for_cellSNP.txt"
    )
    prev_count_file.add_argument('--keep_all_cells', action='store_true',
    help="Use this flag when you want to retain cells not classified"
    " into a donor (e.g. doublets, negatives, etc.)",
    )
    prev_count_file.add_argument('-e', '--extra', nargs='+', 
    help="Classficiation used when a cell is not present "
    "in the classification by hashsolo. NOTE: Accepts multi-word"
    " - space separated input. If not used, use 'None'. "
    "Default: Not Present.", default=["Not", "Present"],
    )
    prev_count_file.add_argument('-d', '--doublet', nargs='?', 
    help="Classficiation used when a cell is a doublet in the"
    " classification by hashsolo. If not used, use 'None'. "
    "Value when parameter present but no value provided: Doublet"
    ". Default: None", const='Doublet', default=None,
    )
    prev_count_file.add_argument('-n', '--negative', nargs='?', 
    help="Classficiation used when a cell is a negative in the "
    "classification by hashsolo. If not used, use 'None'. "
    "Value when parameter present but no value provided: Negative. "
    "Default: None", const='Negative', default=None,
    )

    starsolo_out_file = parser.add_argument_group('mtx_file_path')
    # Optional parameters
    starsolo_out_file.add_argument('--id2name',
    help="File containing 2 columns with gene ids and gene names", 
    )
    starsolo_out_file.add_argument('-m', '--max_mito', nargs='?', type=int,
    help="Max mitochondrial genes(in percent) per cell. Value when "
    "parameter present but no value provided: 5. Default: None", const=5,
    default=None,
    )
    starsolo_out_file.add_argument('-g', '--min_genes', nargs='?', type=int,
    help="Min #genes per cell. Value when parameter present but no value "
    "provided: 1000. Default: None", const=1000, default=None,
    )
    starsolo_out_file.add_argument('--min_cells', nargs='?', type=int, 
    help="Min #cells expressing a gene for it to pass the filter"
    ". Value when parameter present but no value provided: 10"
    ". Default: None", const=10, default=None,
    )
    starsolo_out_file.add_argument('--mito_prefix', nargs='?', help="How mitochondrial"
    " genes can be identified from the gene_info_file. e.g. Value when parameter "
    "present but no value provided'MT-'. Default: None.", const='MT-', default=None,
    )
    
    return parser


def main():
    """Main function"""

    parser = get_argument_parser()
    args = parser.parse_args()

    op=args.output

    # Create parent dir(s) to the output
    if not os.path.isdir(op.replace('/' + os.path.basename(op), '')):
        os.makedirs(op.replace(os.path.basename(op), ''))


    # Previously demultiplexed final count matrix exists
    if args.prev:

        # Join
        extra_val=' '.join(args.extra)

        try:
            adata=ad.read(args.inp)
        except:
            raise ValueError('The cached output(h5ad) of hashsolo doesn\'t exist!')

        # Inform the user when the flag 'keep_all_cells' is used/turned on
        if args.keep_all_cells:
            print("Flag 'keep_all_cells' is used. Hence, all cells will be considered!")
        # Check if the given column exists in the h5ad file
        else:
            if args.column.lower() not in map(str.lower, adata.obs_keys()):
                raise ValueError("Given value for the argument 'column' doesn't exist in the given h5ad file!")

            # No. of Cells
            cells=adata.n_obs
            # Santiy check for the arguments
            # if not args.keep_all_cells:
            if extra_val is not None and extra_val != "None":
                if extra_val.lower() not in adata.obs[args.column].str.lower():
                    raise ValueError("The word(s) provided to the argument 'extra' is not found in the given column (argument 'column')")
                else:
                    extra_val_ser=adata.obs[args.column].str.lower() != extra_val.lower()
            else:
                extra_val_ser=pd.Series(repeat(True, cells))

            if args.doublet.lower() is not None and args.doublet.lower() != "None":
                if args.doublet.lower() not in adata.obs[args.column].str.lower():
                    raise ValueError("The word(s) provided to the argument 'doublet' is not found in the given column (argument 'column')")
                else:
                    doublet_ser=adata.obs[args.column].str.lower() != args.doublet.lower()
            else:
                doublet_ser=pd.Series(repeat(True, cells))

            if args.negative.lower() is not None and args.negative.lower() != "None":
                if args.negative.lower() not in adata.obs[args.column].str.lower():
                    raise ValueError("The word(s) provided to the argument 'negative' is not found in the given column (argument 'column')")
                else:
                    negative_ser=adata.obs[args.column].str.lower() != args.negative.lower()
            else:
                negative_ser=pd.Series(repeat(True, cells))


            adata=adata[ extra_val_ser & doublet_ser & negative_ser ].copy()

    else:
        # Parameters for filtering
        max_mito = args.max_mito
        min_genes = args.min_genes
        min_cells = args.min_cells
        try:
            adata = sc.read_10x_mtx(args.inp[:-13], make_unique=True, var_names= "gene_ids", cache=True)
        except:
            e = sys.exc_info()[0]
            print(f"Error encountered while loading the h5ad file (previously completed demultiplex run)!\nError message: {e}") 

        t2g = pd.read_csv(args.id2name, skiprows=1, usecols=range(2),names=["gene_id", "gene_name"], sep="\t")
        t2g.index = t2g.gene_id
        t2g = t2g.loc[~t2g.index.duplicated(keep='first')]
        adata.var_names_make_unique()
        adata.var["gene_id"] = adata.var.index.values
        adata.var["gene_name"] = adata.var.gene_id.map(t2g["gene_name"])
        adata.var_names = adata.var_names.to_series().map(lambda x: x + '_index')
        adata = adata[:, pd.notna(adata.var["gene_name"])] # Removed gene_ids that don't have an associated gene name
        adata.var_names_make_unique()
        adata.X = adata.X.astype('float64')
        adata.obs_names = adata.obs_names.to_series().map(lambda x: re.sub('-.*', '', x))
        
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        # Filter data wrt mito content
        adata.var["mito"] = adata.var["gene_name"].str.startswith(args.mito_prefix)
        sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
        adata = adata[adata.obs["pct_counts_mito"]< max_mito, :]

    with open(op, 'w') as f:
        for item in adata.obs_names.tolist():
            item=item[:args.barcode_len] # To remove the suffixes used to make barcodes unique, if present
            f.write("%s\n" % item)


if __name__ == '__main__':
    main()