from __future__ import annotations

import argparse
import re

import anndata as ad
import pandas as pd
import scanpy as sc
import scanpy.external as sce


sc.settings.verbosity = 3
sc.logging.print_version_and_date()


def get_argument_parser():
    """Generate and return argument parser."""

    parser = argparse.ArgumentParser(
        description="Run HashSolo demultiplexing on HTO counts"
    )

    parser.add_argument(
        "bustools_out",
        help="Path to HTO h5ad generated from bustools/kite",
    )

    parser.add_argument(
        "matrix_file",
        help="Path to gene expression matrix.mtx.gz",
    )

    parser.add_argument(
        "output_file",
        help="Output h5ad file",
    )

    parser.add_argument(
        "gene_info_file",
        help="Tab-separated gene annotation file",
    )

    # QC parameters
    parser.add_argument(
        "-m",
        "--max_mito",
        type=float,
        default=5,
        help="Maximum mitochondrial percentage",
    )

    parser.add_argument(
        "-g",
        "--min_genes",
        type=int,
        default=1000,
        help="Minimum genes per cell",
    )

    parser.add_argument(
        "-c",
        "--min_cells",
        type=int,
        default=10,
        help="Minimum cells expressing gene",
    )

    parser.add_argument(
        "--mito_prefix",
        default="MT-",
        help="Prefix for mitochondrial genes",
    )

    return parser


def load_gex_data(args):

    t2g = pd.read_csv(
        args.gene_info_file,
        skiprows=1,
        usecols=range(2),
        names=["gene_id", "gene_name"],
        sep="\t",
    )

    t2g.index = t2g.gene_id
    t2g = t2g.loc[~t2g.index.duplicated(keep="first")]

    adata = sc.read_10x_mtx(
        args.matrix_file[:-13],
        make_unique=True,
        var_names="gene_ids",
        cache=True,
    )

    adata.var_names_make_unique()

    adata.var["gene_id"] = adata.var.index.values
    adata.var["gene_name"] = adata.var.gene_id.map(t2g["gene_name"])

    adata.var_names = (
        adata.var_names.to_series().map(lambda x: f"{x}_index")
    )

    adata = adata[:, pd.notna(adata.var["gene_name"])]

    adata.var_names_make_unique()

    adata.obs_names = (
        adata.obs_names.to_series().map(
            lambda x: re.sub(r"-.*", "", x)
        )
    )

    return adata


def perform_qc(
    adata,
    min_genes,
    min_cells,
    max_mito,
    mito_prefix,
):

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata.var["mito"] = adata.var["gene_name"].str.startswith(
        mito_prefix
    )

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mito"],
        inplace=True,
    )

    adata = adata[
        adata.obs["pct_counts_mito"] < max_mito
    ].copy()

    return adata


def main():

    parser = get_argument_parser()
    args = parser.parse_args()

    # Load GEX matrix
    gene_adata = load_gex_data(args)

    # QC
    gene_adata = perform_qc(
        gene_adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mito=args.max_mito,
        mito_prefix=args.mito_prefix,
    )

    # Load HTO counts
    hash_data = ad.read_h5ad(args.bustools_out)

    # Standardize barcode names
    hash_data.obs_names = (
        hash_data.obs_names.to_series().map(
            lambda x: re.sub(r"-.*", "", x)
        )
    )

    # Keep only shared cells
    shared_barcodes = gene_adata.obs_names.intersection(
        hash_data.obs_names
    )

    hash_data = hash_data[shared_barcodes].copy()
    hto_columns = hash_data.obs.columns.tolist()

    sce.pp.hashsolo(hash_data, cell_hashing_columns=hto_columns)

    # Save
    hash_data.write_h5ad(args.output_file)


if __name__ == "__main__":
    main()
    