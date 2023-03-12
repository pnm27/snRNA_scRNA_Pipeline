#!/usr/bin/env python3

# Solo didn't run through scvi, scvi-tools nor scanpy.external
# Only this seems to work
from typing import Union # Need verion > 3.5
import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import os, sys, argparse
from collections import Counter
from collections import defaultdict, OrderedDict as ord_dict
import datetime
from time import sleep
from demultiplex_helper_funcs import demux_by_calico_solo, demux_by_vireo, auto_read


assert sys.version_info >= (3, 5), "This script needs python version >= 3.5!"

# def redo_test(inp_key, rem_last_n_chars=0) -> Union[str, None]:
#     if inp_key in snakemake.input.keys() and os.path.isfile(snakemake.input[inp_key]) and rem_last_n_chars==0:
#         return snakemake.input[inp_key]
#     elif inp_key in snakemake.input.keys() and os.path.isfile(snakemake.input[inp_key]) and rem_last_n_chars>0:
#         return snakemake.input[inp_key][:-rem_last_n_chars]
#     elif inp_key in snakemake.input.keys() and not os.path.isfile(snakemake.input[inp_key]):
#         raise ValueError(f"The file {snakemake.input[inp_key]} does not exist!")
#     else:
#         return None



def get_argument_parser():
    """Generate and return argument parser."""

    #Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Demultiplex pools "
    "(supports hashsolo and vireo)"
    )
    parser.add_argument('input_file', help="Path to matrix.mtx.gz or h5ad "
    "file. If an h5ad file is provided then it is expected that it "
    "has been already processed i.e. poor cells are already filtered out."
    )
    parser.add_argument('count_matrix', help="Path to store the final "
    "count matrix(h5ad)"
    )
    parser.add_argument('gene_info_file', help="Path to file that "
    "contains gene names and ids for annotation (tab-separated txt file)"
    )
    parser.add_argument('--demux_info', help="Path to store demultiplexing "
    "info (tab-separated txt file)"
    )

    # Input of mtx file
    add_redo_grp = parser.add_argument_group('START AFRESH', "Creating output "
    "for the first time."
    )
    add_redo_grp.add_argument('-m', '--max_mito', type=int, help="Max "
    "mitochondrial genes(in percent) per cell. Default: 5", 
    default=5
    )
    add_redo_grp.add_argument('--mito_prefix', help="Prefix for mitochondrial "
    "genes. Default: 'MT-'", default="MT-" 
    )
    add_redo_grp.add_argument('-g', '--min_genes', type=int, help="Min #genes "
    "per cell. Default: 1000", default=1000
    )
    add_redo_grp.add_argument('-c', '--min_cells', type=int, help="Min #cells "
    "expressing a gene for it to pass the filter. Default: 10", default=10
    )  

    # For calico_solo inputs
    cs = parser.add_argument_group('HASHSOLO DEMUX OPTIONS', "Add calico_solo "
    " demultiplex to create final count matrix file")
    cs.add_argument('-w', '--wet_lab_file', help="Path to file that "
    "contains HTO info for each set (either csv or tsv file)"
    )
    cs.add_argument('--calico_solo', dest='hashsolo_out', help="Path "
    "to cached output of hashsolo(h5ad)", metavar="hashsolo.h5ad"
    )
    cs.add_argument('--hto_sep', help="If, per each sample in the "
    "wet lab file (6th positional argument to this script), HTOs are "
    "all present in one row separated by some SEP then specify it here. "
    "Default: ' '", default=' '
    )
    cs.add_argument('--columns', nargs=4, help="List of column names "
    "RESPECTIVELY to sample_ID (as present in the spreadsheets - "
    "identifies pooled samples), HTO numbers, it's associated barcodes, "
    "and Donors/SubIDs (contains each multiplexed donors).", 
    metavar=('sample_ID', 'HTO_name', 'HTO_barcode', 'Sub_ID'),
    default=['unique_sample_ID', 'hashtag', 'ab_barcode', 'SubID']
    )
    cs.add_argument('-s', '--sample_name', help="Name of the sample. "
    "Check whether this name coincides with that in this script as well "
    "as the one in the wet_lab_file"
    )
    cs.add_argument('--no-demux-stats-cs', action='store_true', 
            dest="cs_stats",
			help="If flag is used no demux stats are present",
			)
    cs.add_argument('--no-subid_convert', action='store_true', 
            dest="subid_convert",
			help="If flag is used no conversion to subID is needed",
			)

    # For vireo inputs
    vs = parser.add_argument_group("VIREO DEMUX OPTIONS", "Add "
    "genotype-based demultiplexing outputs to create final count matrix"
    )
    vs.add_argument('--vireo_out', help="Path to donor_ids.tsv file",
    metavar="donor_ids.tsv")
    vs.add_argument('--converter_file', help="If names from vireo output "
    "needs to be changed"
    )
    vs.add_argument('--no-demux-stats-vs', action='store_true',
            dest="vs_stats",
			help="If flag is used no demux stats are present",
			)

    return parser


def main():
    """Main entry point"""

    parser = get_argument_parser()
    args = parser.parse_args()

    sc.settings.set_figure_params(dpi_save=400, format='png', 
                                color_map = 'viridis_r')
    sc.settings.autosave = True
    sc.settings.autoshow = False
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_version_and_date()

    # DEPRACATED-------------------------------------------------------------
    # Adding demultiplex info or creating "new" final count matrix 
    # redo = True if args.matrix_file is None else False

    # If 'redoing' demultiplexing then whether vireo output is to be added 
    # or calico_solo output is to be added to an existing final count matrix file
    # add_calico = redo_test('vireo_out') if redo is None else None
    # add_vireo = redo_test('calico_solo_out') if redo is None else None


    # If creating "new" final count matrix then selecting both or one
    # starsolo_mat = redo_test('starsolo_out', 13) if redo is None else None #args.matrix_file[:-13]
    # calico_demux = redo_test('calico_solo_out') if redo is None else None
    # vireo_demux = redo_test('vireoSNP_out') if redo is None else None

    # -----------------------------------------------------------------------
    # If an mtx file is provided then it is expected to be not-filtered
    # i.e. bad quality cells and genes exist
    # While an h5ad file is assumed to be rid of these issues
    redo = True if args.input_file.endswith('.h5ad') else False

    add_calico = args.hashsolo_out
    add_vireo = args.vireo_out

    starsolo_mat = args.input_file[:-13] if not redo else None

    # Filtering stats
    filter_info = []
    # Individual demux stats
    solo_dem_stats = []
    vir_dem_stats = []

    # Both calico_solo and vireo can't be None
    if add_calico is None and add_vireo is None:
        raise ValueError(
            "For this script provide either a VALID calico_solo h5ad "
            "or vireoSNP's donor file"
            )
    
    no_stats = True if redo and args.cs_stats and args.vs_stats else False
    
    # For assigning gene names
    t2g = pd.read_csv(args.gene_info_file, skiprows=1, usecols=range(2),
                    names=["gene_id", "gene_name"], sep="\t")
    t2g.index = t2g.gene_id
    t2g = t2g.loc[~t2g.index.duplicated(keep='first')]
    

    # Store output_file and Create necessary folders
    op = args.count_matrix

    # Create parent dir(s) to the output
    if not os.path.isdir(op.replace('/' + os.path.basename(op), '')):
        os.makedirs(op.replace(os.path.basename(op), ''))

    # Initial run of demultiplexing
    if not redo:
        # batch for wet lab file 
        # samp=batch

        # Parameters for filtering
        max_mito=args.max_mito
        min_genes=args.min_genes
        min_cells=args.min_cells

        # Process STARsolo output----------------------------------------------------------------------------------------------------------------------------
        ct = datetime.datetime.now()
        print(f"Processing STARsolo's (or mtx files) output at: {ct}")
        try:
            adata = sc.read_10x_mtx(starsolo_mat, make_unique=True, 
            var_names= "gene_ids", cache=True)
        except:
            e = sys.exc_info()[0]
            print("Error encountered while loading the mtx files!\nError "
                f"message: {e}")

        
        print(adata)
        filter_info.append(( 'Started with cells', adata.n_obs))
        filter_info.append(( 'Started with genes', adata.n_vars))
        adata.var_names_make_unique()
        adata.var["gene_id"] = adata.var.index.values
        adata.var["gene_name"] = adata.var.gene_id.map(t2g["gene_name"])
        adata.var_names = (
            adata.var_names.to_series().map(lambda x: x + '_index')
            )
        total_umi_lost = adata[:, ~ pd.notna(adata.var["gene_name"])].X.sum()
        avg_umis_lost_per_cell = (
            adata[:, ~ pd.notna(adata.var["gene_name"])].X.sum(axis=1)
            .mean()
            )
        avg_umis_lost_per_gene = (
            adata[:, ~ pd.notna(adata.var["gene_name"])].X.sum(axis=0)
            .mean()
            )
        # Removed gene_ids that don't have an associated gene name
        adata = adata[:, pd.notna(adata.var["gene_name"])]
        adata.var_names_make_unique()
        adata.X = adata.X.astype('float64')
        filter_info.append(( 'gene_ids with an associated gene_name', 
                            adata.n_vars))
        filter_info.append(( 'total UMI counts lost to gene ids wo names', 
                            total_umi_lost))
        filter_info.append(( 'avg UMI counts lost per cell', 
                            avg_umis_lost_per_cell))
        filter_info.append(( 'avg UMI counts lost per gene', 
                            avg_umis_lost_per_gene))
        print(adata)

        # Filter data using cell level metrics
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        filter_info.append(( 'min #genes expressed per cell', 
                            min_genes))
        filter_info.append(( 'min #cells expressing per gene', 
                            min_cells))
        filter_info.append(( 'Retained cells after previous filter', 
                            adata.n_obs))
        filter_info.append(( 'Retained genes after previous filter', 
                            adata.n_vars))

        # Filter data wrt mito content
        adata.var["mito"] = adata.var["gene_name"].str.startswith(args.mito_prefix)
        sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
        adata = adata[adata.obs["pct_counts_mito"]< max_mito, :]
        filter_info.append(( 'max percent mito content per cell', 
                            max_mito))
        filter_info.append(( 'cells with low mito percent', 
                            adata.n_obs))
        print(adata)
    
    else:
        ct = datetime.datetime.now()
        print(f"Reading given h5ad input file at: {ct}")
        try:
            adata = ad.read(args.input_file)
        except:
            e = sys.exc_info()[0]
            print(
                "Error encountered while loading the input h5ad file!"
                f"\nError message: {e}"
                )

        print("Successfully loaded the input file!")



    # Batch info
    # This is the values that will be stored in the final h5ad file
    # batch=args.sample_name.replace('-', '_')+'_cDNA'
    # Prepare for Extra Information
    replicate=args.sample_name.split('_')[2]
    # add few more annotations
    adata.obs['batch'] = args.sample_name
    adata.obs['rep'] = replicate
    adata.obs['set'] = '_'.join(args.sample_name.split('_')[:3])[:-1]

    # For demultiplexing using calico_solo
    if add_calico is not None:

        ct = datetime.datetime.now()
        print(
            f"Starting calico_solo/hashsolo demultiplexing at: {ct}"
        )

        cols = args.columns

        # Wet Lab file, Filter wet lab file's columns, if needed
        df = auto_read(args.wet_lab_file)
        if df.loc[df[cols[0]].str.lower() == args.sample_name.lower()].empty:
            raise ValueError("Check dtypes!\nSample (variable name 'var'"
            f", data type {type(args.sample_name)}, with value "
            f"{args.sample_name.lower()} ) couldn't be subset from the wet lab "
            "file.\nData types for the wet lab "
            f"file:\n{df.dtypes}\nWhile the top 5 rows are:\n{df[cols[0]][:5]}")
        else:
            df = df.loc[df[cols[0]].str.lower() == args.sample_name.lower()]


        # Load hashsolo/calico_solo output (h5ad)
        try:
            dem_cs = ad.read(add_calico)
        except:
            e = sys.exc_info()[0]
            print(
                "Error encountered while loading the h5ad file!"
                f"\nError message: {e}"
                )

        print("Successfully loaded calico_solo output!")

        # Demultiplex and assign samples-------------------------------------

        # hto_tags_cs, hto_tags_ms, and hto_tags_hd contain (in sequence): 
        # pd DF with barcodes as index. SubID and HTO number (HTO1, HTO2, etc)
        #  as columns no. of doublet cells, no. of negative cells]

        ct = datetime.datetime.now()
        print(
            "Starting: Assigning cell classifications by"
            f"hashsolo/calico solo at: {ct}"
            )
        cs_dons, hto_name_cs, temp_df = demux_by_calico_solo(
            adata.obs_names.to_series(), df, args.sample_name, 
            args.hto_sep, [cols[1], cols[3]], 
            dem_cs.obs['Classification'], args.subid_convert
            )
        
        ct = datetime.datetime.now()
        print(
            "Assigning demultiplexing info for hashsolo/calico "
            f"solo at: {ct}"
            )

        adata.obs['SubID_cs'] = cs_dons
        adata.obs['HTO_n_cs'] = hto_name_cs
        if not args.cs_stats:
            solo_dem_stats.extend(temp_df)

        # If Subject IDs aren't 'string' then convert them
        adata.obs['SubID_cs']=adata.obs['SubID_cs'].apply(str)

    # To do
    if add_vireo is not None:
        print("Starting demultiplexing through vireoSNP's output")
        # Demultiplex and assign samples------------------------------------

        ct = datetime.datetime.now()
        print(
            "Starting: Assigning cell classifications by"
            f"vireoSNP at: {ct}"
            )
        vs_dons, temp_df = demux_by_vireo(
            adata.obs_names.to_series(), add_vireo, args.converter_file
            )
        
        ct = datetime.datetime.now()
        print(
            "Assigning demultiplexing info for vireo "
            f"solo at: {ct}"
            )

        adata.obs['SubID_vs'] = vs_dons

        if not args.vs_stats:
            vir_dem_stats.extend(temp_df)
        
        
 

    ct = datetime.datetime.now()
    print(
        "Saving All demultiplex info as a tsv file and also the h5ad files: "
        f"{ct}"
        )
    
    if args.demux_info is not None:
        # Need only cs stats
        if solo_dem_stats and not vir_dem_stats:
            filter_info.extend(solo_dem_stats)
        elif not solo_dem_stats and vir_dem_stats:
            filter_info.extend(vir_dem_stats)
        else:
            print("This error shouldn't occur!")

        solo_run_df = pd.DataFrame(filter_info, 
                                    columns=['Observations', 'Vals'])
        
        solo_run_df.to_csv(args.demux_info, sep = "\t", index=False)

    adata.write(op)

    ct = datetime.datetime.now()
    print(f"Finished: Processing Sample {batch} at: {ct}")



if __name__ == '__main__':

    main()
    sleep(60)