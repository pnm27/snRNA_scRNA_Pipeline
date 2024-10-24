#!/usr/bin/env python3
import pandas as pd, anndata as ad
import os, re, argparse
from collections import Counter


# Save dataframe to file
def save_df(df, suff, op, overwrite=False):
    
    if overwrite and os.path.isfile(op):
        os.remove(op)
            
    if suff.lower() is None:
        df.to_csv(op, sep = "\t", index=False)
    elif suff.lower() == '.csv':
        df.to_csv(op, index=False)
    elif suff.lower() == '.tsv' or suff.lower() == '.txt':
        df.to_csv(op, sep = "\t", index=False)
    else:
        print("Can't save to file because of incorrect extension. "
    "Extension can be 'csv', 'tsv' or 'txt'")


def get_argument_parser():
    """Generate and return argument parser."""

    #Parse Command-Line arguments
    parser = argparse.ArgumentParser(
        description="Create a two-columned txt file with barcodes "
        "corresponding to each Donor when provided with annotated h5ad "
        "or donor_ids.tsv (vireo output)."
        )
    parser.add_argument('inp', help="Path to cached final count matrix "
            "with annotated barcodes or donor_ids.tsv (vireo output)")
    parser.add_argument('output', help="Provide output file name with "
            "extension. Supports only csv, tsv or txt (tab sep).")
    
    # Optional parameters

    parser.add_argument('--overwrite', action='store_true', 
            help="Flag describing whether to overwrite or not. "
            "If not, \"_2\" will be appended."
            )
    parser.add_argument('--converter_file', 
            help="If names in the h5ad classification or donor_ids.tsv "
            " need to be changed. Only csv file is supported."
            )
    parser.add_argument('--converter_file_from_col', 
            help="Column containing matching, with respect to that in "
            "vireo output or h5ad donor columns, donor names.", dest='from_col',
            )
    parser.add_argument('--converter_file_to_col', 
            help="Column containing the final donor names", dest='to_col',
            )
    
    # subparsers = parser.add_subparsers(help='subcommand help', required=True, 
    #         dest='subparser_name')
    
    # create the parser for the "demux_h5ad" command
    # demux_h5ad = subparsers.add_parser('demux_h5ad', help="Demux h5ad using"
    #         " columns containing classifications."
    #         )
    
    parser.add_argument('--h5ad_donor_column', 
            help="If the input is an h5ad containing donor classification, "
            "provide the 'obs' column containing the classification", 
            dest='h5ad_don_col',
            )
        
    parser.add_argument('--split_by', nargs='+', 
            help="To split bams by classifications preoduced by either "
            "one or both of calico_solo and vireo provide column(s) "
            "that contain the classification. If vireo, specify 'vireo'. "
            "Use when a h5ad containing classification is provided as input.",
            dest="method",
            )
    parser.add_argument('--demux_suffix', nargs='*', 
            help="Suffix(es) to use for naming the output files according "
            "to the demultiplexing methods specified.",
            metavar="suffix",
            )
    # demux_h5ad.set_defaults(func=h5ad_demux)
    # USEFUL WHEN SPLITTING BAMS BY MULTIPLE DEMUX METHODS

    
    # create the parser for the "demux_vireo" command
    # demux_vireo = subparsers.add_parser('demux_vireo', help='a help')
    # demux_vireo.set_defaults(func=demuxByVireo)

    return parser


def main():
    r"""Main entry function

    This function creates a 2-columned text file with donors as the 
    first column with its corresponding barcodes in the second.
    NOT SUPPORTED YET:
    If multiple demultiplexing softwares have annotated the cells
    and one wants to create sep output files for each of them
    """

    parser = get_argument_parser()
    # Parse arguments
    args = parser.parse_args()

    inp_h5ad = False if os.path.basename(args.inp) == 'donor_ids.tsv' else True
    fout = args.output
    rem_op = args.overwrite

    

    if args.converter_file is not None:
        conv_df = pd.read_csv(args.converter_file)
        assert args.from_col in conv_df.columns, "The provided 'from' column is not present in the converter file!"
        assert args.to_col in conv_df.columns, "The provided 'to' column is not present in the converter file!"

    # USEFUL FOR MULTI DEMUX METHOD CASE
    file_ext = re.search(r'(\.[^.]+)$', fout)
    file_basename = fout.replace(file_ext.group(1), '') if file_ext is not None else fout

    # Call the respective functions
    # args.func(args)

    # This is for splitting the bam files after finalizing genotypes
    if inp_h5ad:
        adata = ad.read_h5ad(args.inp)
        assert args.h5ad_don_col in adata.obs.columns, ("The provided h5ad 'donor' column is not "
            "present in the h5ad file!")
        # For donor name conversion
        if args.converter_file is not None:
            conv_df = conv_df.loc[ conv_df[args.from_col].isin(adata.obs[args.h5ad_don_col].unique()), 
                                    [args.to_col, args.from_col] ]
            conv_map = conv_df.set_index(args.from_col)[args.to_col].to_dict()

        # Support for calico_solo based demux too
        demux_methods = args.method if isinstance(args.method, list) else [args.method]
        # # USEFUL WHEN SPLITTING BAMS BY MULTIPLE DEMUX METHODS
        if len(demux_methods) > 1:
            assert len(demux_methods) == len(args.demux_suffix), "Number of suffixes"
            " and number of demux methods aren't the same!"
        suff = args.demux_suffix[1:] if args.demux_suffix.startswith('_') else args.demux_suffix

        for i, dem in enumerate(demux_methods):

            d = Counter(adata.obs[args.h5ad_don_col])

            if dem in ["vs", "vireo"]:
                for vals in ['Doublet', 'Negative', 'Not Present']:
                    d.pop(vals, None)
            elif dem in ["cs", "calico", "calico_solo", "hashsolo"]:
                for vals in ['doublet', 'unassigned']:
                    d.pop(vals, None)


            df_l = []
            for k in d.keys():
                df_l.extend( (k, bc) for bc in adata[adata.obs[args.h5ad_don_col] == k].obs_names.to_list() )

            temp_df = pd.DataFrame(df_l, columns=['Subj_ID', 'barcodes'])
            # For donor name conversion
            if args.converter_file is not None:
                temp_df.loc[:, 'Subj_ID'] = temp_df.loc[:, 'Subj_ID'].astype('category')
                temp_df.loc[:, 'Subj_ID'] = temp_df.loc[:, 'Subj_ID'].map(conv_map)
            
            fout = file_basename + suff[i] + file_ext if len(demux_methods) > 1 else fout
            save_df(temp_df, file_ext, fout)

    else:
        vir_class = pd.read_csv(args.inp, sep='\t', usecols=["cell", "donor_id"])
        vir_class = vir_class[(vir_class["donor_id"] != "doublet") & (vir_class["donor_id"] != "unassigned")]
        vir_class.reset_index(drop=True, inplace=True)

        # For converting genotype IDs to Subject IDs-------------------------------
        if args.converter_file is not None:
            # Project-specific filtering of donor_ids for conversion
            # e.g.
            # vir_class["donor_id"] = vir_class["donor_id"]
            #                           .apply(lambda x: '_'.join(x.split('_')[1:]) \
            #                                 if x.startswith('0_') else x)

            
            conv_df = conv_df.loc[conv_df[args.from_col].isin(vir_class['donor_id'].unique()), [args.to_col, args.from_col]]
            conv_map = conv_df.set_index(args.from_col)[args.to_col].to_dict()
            # Keep unknown donors the same (unknown donors will be called donor0, donor1, etc.)
            conv_map.update({ k:k for k in vir_class['donor_id'].unique() if k.startswith('donor') })
            
            vir_class['donor_id'] = vir_class['donor_id'].astype('category')
            vir_class['donor_id'] = vir_class['donor_id'].map(conv_map)
            vir_class = vir_class.loc[:, ["donor_id", "cell"]]
            
            vir_class.rename(columns={"cell":"barcodes", "donor_id":"Subj_ID"}, inplace=True, errors="raise")
            
        # -------------------------------------------------------------------------
        else:
            vir_class.rename(columns={"cell":"barcodes", "donor_id":"Subj_ID"}, inplace=True, errors="raise")

        vir_class = vir_class.loc[:, ["Subj_ID", "barcodes"]]

        save_df(vir_class, file_ext, fout)


    # # OLDER CODE ----------------------------
    # if rem_op:
    #     try:
    #         os.remove(op_f)

    #     except:
    #         pass
    #     save_df(temp_df, file_ext, op_f)
    # # ---------------------------------------------
    # else:   
    #     file_pref = file_basename + '_' + suff + '_2'
    #     new_name = file_pref + '_2' + file_ext
    #     save_df(temp_df, file_ext, new_name)
    # # ---------------------------------------------------------------------------


if __name__ == '__main__':
    main()
    