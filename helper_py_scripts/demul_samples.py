#!/usr/bin/env python3

# Solo didn't run through scvi, scvi-tools nor scanpy.external
# Only this seems to work
from typing import Union # Need verion > 3.5
import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import os, sys, argparse, itertools
from collections import Counter
from collections import defaultdict, OrderedDict as ord_dict
import datetime
from time import sleep
from demultiplex_helper_funcs import (
    auto_read, 
    demux_by_calico_solo, 
    demux_by_vireo,
    get_donor_info, 
    )


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

# Example run: For multi-HTO per donor in a multiplexed pool

# python3 demul_samples.py [--demux_info demultiplex/info/poolA.txt -m 5 --mito_prefix 5 -g 1000 -c 10 -w wet_lab_file.csv
                        # --calico_solo demultiplex/solo/poolA/round1.h5ad demultiplex/solo/poolA/round2.h5ad 
                        # --hto_sep "_" --columns Pool_name HTO HTO_bc Sample -p poolA
                        # -pref round1 round2 [--no-demux-stats-cs] [--no-subid_convert] 
                        # input_file count_matrix gene_info_file

# 
# wet_lab_file.csv
# Pool_name Sample HTO HTO_bc
# poolA	samp1	HTO1_HTO2	AATCGATCGAT_ATCGGCTTAGCT
# poolA	samp2	HTO7_HTO8	AATCGGCTGAT_AAACGCTTAGCT
# poolA	samp3	HTO3_HTO5	AACCCATCGAT_ATCGCTCCAGCT
# poolA	samp4	HTO3_HTO4	AAGCGATCGAT_ATTTCGCTTAGC
# poolA	samp5	HTO5_HTO6	AATCGACTGAT_ATCGGCTTACCC
# poolA	samp6	HTO9_HTO10	AATCCTCCGAT_ATCGGCTGGGC
# poolB	samp7	HTO7_HTO8	AATCGGCTGAT_AAACGCTTAGCT
# poolB	samp8	HTO9_HTO10	AATCCTCCGAT_ATCGGCTGGGC
# poolB	samp9	HTO6_HTO5	AATCGACTGAT_ATCGGCTTACCC
# poolB	samp10	HTO2_HTO5	ATCGGCTTAGCT_ATCGGCTTACCC
# poolB	samp11	HTO1_HTO2	AATCGATCGAT_ATCGGCTTAGCT
# poolB	samp12	HTO3_HTO4	AAGCGATCGAT_ATTTCGCTTAGC

# and calico_solo output 'tree' looks like:

# demultiplex/
# └── solo
#     ├── poolA
#     │   ├── round1.h5ad
#     │   └── round2.h5ad
#     └── poolB
#         ├── round1.h5ad
#         └── round2.h5ad




def get_argument_parser():
    """Generate and return argument parser."""

    #Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Demultiplex pools "
    "(supports hashsolo and vireo). Note: poor cells are not removed but "
    "they aren't included while demultiplexing"
    )
    parser.add_argument('input_file', help="Path to matrix.mtx.gz or h5ad "
    "file. If an h5ad file is provided then it is expected that it "
    "has been already processed i.e. poor cells are already filtered out."
    )
    parser.add_argument('count_matrix', help="Path to store the final "
    "count matrix(h5ad)"
    )
    parser.add_argument('gene_info_file', help="Path to the file that "
    "contains gene names and ids for annotation (tab-separated txt file)"
    )
    parser.add_argument('--demux_info', help="Path to store demultiplexing "
    "info (tab-separated txt file)"
    )

    # Input of mtx file
    add_redo_grp = parser.add_argument_group('START AFRESH', "Creating output "
    "for the first time."
    )
    add_redo_grp.add_argument('-m', '--max_mito', type=int, nargs='?', 
    help="Max mitochondrial genes(in percent) per cell. "
    "If no given value to parameter, will default to 5 otherwise no filter", 
    default=None, const=5,
    )
    add_redo_grp.add_argument('--mito_prefix', nargs='?', 
    help="Prefix for mitochondrial genes. If no given value to "
    "parameter, will default to 'MT-' otherwise no filter", 
    default=None, const="MT-",
    )
    add_redo_grp.add_argument('-g', '--min_genes', nargs='?', type=int, 
    help="Min #genes per cell. If no given value to parameter, will "
    "default to 1000 otherwise no filter", default=None, const=1000,
    )
    add_redo_grp.add_argument('-c', '--min_cells', nargs='?', type=int, 
    help="Min #cells expressing a gene for it to pass the filter. "
    "If no given value to parameter, will default to 10 otherwise "
    "no filter", default=None, const=10,
    )  

    # For calico_solo inputs
    cs = parser.add_argument_group('HASHSOLO DEMUX OPTIONS', "Add calico_solo "
    " demultiplex to create final count matrix file")
    cs.add_argument('-w', '--wet_lab_file', help="Path to file that "
    "contains HTO info for each set (either csv or tsv file)"
    )
    cs.add_argument('--calico_solo', dest='hashsolo_out', help="Path "
    "to cached output of hashsolo(h5ad) file(s). If no given value to "
    "parameter, will default to not process hashsolo output otherwise ",
    nargs='*', metavar="hashsolo.h5ad", default=None,
    )
    cs.add_argument('--hto_sep', nargs='?', help="If, per each pool in the "
    "wet lab file (6th positional argument to this script), HTOs are "
    "all present in one row separated by some SEP then specify it here. "
    "NOTE: A value is also expected"
    "Default: None", default="_", const=None,
    )
    cs.add_argument('--columns', nargs=4, help="List of column names "
    "RESPECTIVELY to sample_ID (as present in the spreadsheets - "
    "identifies pooled samples), HTO numbers, it's associated barcodes, "
    "and Donors/SubIDs (contains each multiplexed donors). If multiple HTOs "
    "are used per samples in a pool then use 'hto_sep' parameter "
    "to explain how the htos are present. Check 'hto_sep' description "
    "for more info", 
    metavar=('pool_ID', 'HTO_name', 'HTO_barcode', 'donor_ID'),
    default=['unique_sample_ID', 'hashtag', 'ab_barcode', 'SubID']
    )
    cs.add_argument('-p', '--pool_name', help="Name of the pool. "
    "Check whether this name coincides with that in this script as well "
    "as the one in the wet_lab_file"
    )
    cs.add_argument('--no-demux-stats-cs', action='store_true', 
            dest="cs_stats",
			help="If flag is used no demux stats will be stored.",
			)
    cs.add_argument('--no-subid_convert', action='store_true', 
            dest="subid_convert",
			help="If flag is used no conversion to subID is needed."
            "Also expected when used for multi-HTO setup",
			)
    # Only for multi-HTO pools
    cs.add_argument('--pref', help="Prefix for each hashsolo "
    "run per pool. This should match the number of hashsolo files "
    "provided as input and in the same sequence! ", nargs='*', metavar="prefix",
    )

    # For vireo inputs
    vs = parser.add_argument_group("VIREO DEMUX OPTIONS", "Add "
    "genotype-based demultiplexing outputs to create final count matrix."
    )
    vs.add_argument('--vireo_out', help="Path to donor_ids.tsv file",
    metavar="donor_ids.tsv"
    )
    vs.add_argument('--converter_file', help="If names from vireo output "
    "needs to be changed."
    )
    vs.add_argument('--no-demux-stats-vs', action='store_true',
            dest="vs_stats",
			help="If flag is used no demux stats will be stored.",
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

    multi_hto_setup = True if len(add_calico) > 1 else False

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
    
    # no_stats = True if args.cs_stats and args.vs_stats else False
    
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
        # Typically shows discrepancies in reference genome used for
        # alignemnt and that for annotation (here: gene_info file)
        geneids_w_gene_names = pd.notna(adata.var["gene_name"]) 
        adata.var_names = (
            adata.var_names.to_series().map(lambda x: x + '_index')
            )
        total_umi_lost = adata[:, ~ geneids_w_gene_names].X.sum()
        # avg_umis_per_cell_before = (
        #     adata.X.sum(axis=1)
        #     .mean()
        #     )
        # avg_umis_per_gene_before = (
        #     adata.X.sum(axis=0)
        #     .mean()
        #     )
        avg_umis_lost_per_cell = (
            adata[:, ~ geneids_w_gene_names].X.sum(axis=1)
            .mean()
            )
        avg_umis_lost_per_gene = (
            adata[:, ~ geneids_w_gene_names].X.sum(axis=0)
            .mean()
            )
        
        # Don't consider gene_ids that don't have an associated gene name      
        adata.var_names_make_unique()
        adata.X = adata.X.astype('float64')
        filter_info.append(( 'gene_ids with an associated gene_name', 
                            adata.n_vars))
        filter_info.append(( 'avg UMI counts lost per cell', 
                            avg_umis_lost_per_cell))
        filter_info.append(( 'avg UMI counts lost per gene', 
                            avg_umis_lost_per_gene))
        filter_info.append(( 'total UMI counts lost to gene ids wo names', 
                            total_umi_lost))
        print(adata)

        # Filter data using cell level metrics
        cell_subset, n_genesPerCell = sc.pp.filter_cells(adata, min_genes=min_genes, inplace=False)
        adata.obs['n_genes'] = n_genesPerCell
        gene_subset, n_cellsPerGene = sc.pp.filter_genes(adata, min_cells=min_cells, inplace=False)
        adata.var['n_cells'] = n_cellsPerGene
        filter_info.append(( 'min #genes expressed per cell', 
                            min_genes))
        filter_info.append(( 'min #cells expressing per gene', 
                            min_cells))
        filter_info.append(( 'Remaining cells after previous filter', 
                            cell_subset.sum()))
        filter_info.append(( 'Remaining genes after previous filter', 
                            gene_subset.sum()))

        # Filter data wrt mito content
        adata.var["mito"] = adata.var["gene_name"].str.startswith(args.mito_prefix)
        sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
        mito_QCpass = adata.obs["pct_counts_mito"]< max_mito
        # adata._inplace_subset_obs(cell_subset)
        # adata._inplace_subset_obs(mito_QCpass)
        # adata._inplace_subset_var(gene_subset)
        filter_info.append(( 'max percent mito content per cell', 
                            max_mito))
        filter_info.append(( 'cells with low mito percent', 
                            mito_QCpass.sum()))
        
        qc_pass_cells = cell_subset.astype('bool') & mito_QCpass.astype('bool')
        filter_info.append(( 'Cells passing mito and basic filter threshold', 
                            qc_pass_cells.sum()))
        
        # Assign 'QC_pass' observation according to the selected filters
        adata.obs['QC_pass'] = qc_pass_cells
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
    # batch=args.pool_name.replace('-', '_')+'_cDNA'
    # Prepare for Extra Information
    replicate=args.pool_name.split('_')[2]
    # add few more annotations
    adata.obs['batch'] = args.pool_name
    adata.obs['rep'] = replicate
    adata.obs['set'] = '_'.join(args.pool_name.split('_')[:3])[:-1]
    cell_bcs =  adata.obs_names.to_series()

    # For demultiplexing using calico_solo
    if add_calico is not None:
        # Wet Lab file, Filter wet lab file's columns, if needed
        df = auto_read(args.wet_lab_file)
        if df.loc[df[cols[0]].str.lower() == args.pool_name.lower()].empty:
            raise ValueError("Check dtypes!\nSample (variable name 'var'"
            f", data type {type(args.pool_name)}, with value "
            f"{args.pool_name.lower()} ) couldn't be subset from the wet lab "
            "file.\nData types for the wet lab "
            f"file:\n{df.dtypes}\nWhile the top 5 rows are:\n{df[cols[0]][:5]}")
        else:
            df = df.loc[df[cols[0]].str.lower() == args.pool_name.lower()]

        ct = datetime.datetime.now()
        print(
                f"Starting calico_solo/hashsolo demultiplexing at: {ct}"
            )

        cols = args.columns
        hto_count=0

        # If more than 1 hasholo outputs are given then it means
        # more than 1 HTO is used for sample identification
        if multi_hto_setup:
            assert args.hto_sep is not None, "Expected command line " \
            "argument to 'hto_sep' parameter as multiple calico_solo " \
            "outputs are provided per pool!"
            assert args.subid_convert == False, "'no-subid-convert' flag " \
            "is not expected as multiple calico_solo outputs are " \
            "provided per pool!"
            assert len(add_calico) == len(args.pref), "Prefix options " \
            "given to 'pref' parameter should be equal to number of " \
            "hashsolo inputs!"

        for c in add_calico:
            # For use in column names
            suff = args.pref

            print(
                f"Loading file {c}"
            )
        
            # Load hashsolo/calico_solo output (h5ad)
            try:
                dem_cs = ad.read(c)
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
            if not multi_hto_setup:
                cs_dons, hto_name_cs, temp_df = demux_by_calico_solo(
                    cell_bcs, df, args.pool_name, 
                    args.hto_sep, [cols[1], cols[3]], 
                    dem_cs.obs['Classification'], args.subid_convert,
                    # hto_count, multi_hto_setup
                    )
            else:
                cs_dons, hto_name_cs, temp_df = demux_by_calico_solo(
                    cell_bcs, df, args.pool_name, 
                    args.hto_sep, [cols[1], cols[3]], 
                    dem_cs.obs['Classification'], False,
                    # hto_count, multi_hto_setup
                    )
            
            ct = datetime.datetime.now()
            print(
                "Assigning demultiplexing info for hashsolo/calico "
                f"solo at: {ct}"
                )

            
            # When not a multi-hto experiment
            if not multi_hto_setup and not args.subid_convert:
                adata.obs['SubID_cs'] = cs_dons
                adata.obs['HTO_n_cs'] = hto_name_cs
            #  No donor-name conversion needed
            elif not multi_hto_setup and args.subid_convert:
                adata.obs['HTO_n_cs'] = hto_name_cs
            # For multi-hto setup donor-names not needed per each run 
            # of calico solo
            else:
                adata.obs['HTO_n_cs_'+suff] = hto_name_cs

            if not args.cs_stats and not multi_hto_setup:
                solo_dem_stats.extend(temp_df)

            

        # args.subid_convert is assumed to be False i.e. in multi-HTO setup
        # hto combos need to be changed to donor names
        if multi_hto_setup:
            hto_cols = [c for c in adata.obs.columns if c.startswith('HTO_n_cs_')]
            adata.obs.loc[:, 'combo'] = adata.obs[hto_cols[0]].str.cat(
                [adata.obs[x] for x in hto_cols[1:]], sep='_'
                )
            
            temp_l = []
            for i, j in enumerate(df[cols[1]]):
                vals = [df[cols[3]][i], j]
                vals.extend([list(c) \
                    for c in list(set(list(itertools.permutations(
                    j.split(args.hto_sep)))))
                    ]
                )
                temp_l.append(vals)
            
            temp_df2 = pd.DataFrame(temp_l, 
                        columns=['comb' + str(x+1) for x in range(len(vals)-2)]
                        )

            adata.obs['SubID_cs'] = adata.obs['comb'].apply(get_donor_info, args=(temp_df2, ))

        # If Subject IDs aren't 'string' then convert them
        if 'SubID_cs' in adata.obs.columns:
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
            f"at: {ct}"
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
        elif solo_dem_stats and vir_dem_stats:
            filter_info.extend(solo_dem_stats)
            filter_info.extend(vir_dem_stats)
        else:
            print("Proivded a file to collect demultiplex info but used "
                  "flags to not include both demux stats!")

        solo_run_df = pd.DataFrame(filter_info, 
                                    columns=['Observations', 'Vals'])
        
        solo_run_df.to_csv(args.demux_info, sep = "\t", index=False)

    adata.write(op)

    ct = datetime.datetime.now()
    print(f"Finished: Processing Sample {args.pool_name} at: {ct}")



if __name__ == '__main__':

    main()
    sleep(60)