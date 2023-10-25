#!/usr/bin/env python3

# Example:
# If h5ad files present in the dir "/sc/arion/projects/CommonMind/pnm/AMP_PD/final_count_matrix/test/"
# need to be aggregated with the help of 2 conversion tables, which are as follows:
# 1) conversion file 1: conversion_vir_w_gt2.tsv
#     Sample Donor_name final_vcfID donor Final_breg
#     PD-set1-1 246041 246041 246044 BM-17
#     PD-set1-1 246048 246048 243003 BM-4 
#     PD-set1-1 donor2 250118  
#     PD-set1-2 246041 246041 246044 BM-17
#     PD-set1-2 donor1 250118  
#     PD-set1-2 246048 246048 243003 BM-4
#     ...

# 2) conversion file 2: conversion_vir_irst_run2.tsv
#     Sample Donor_name donID donor
#     PD-set1-1 donor0 250118 250149 
#     PD-set1-1 donor1 246048 243003 
#     PD-set1-1 donor2 246041 246044 
#     PD-set1-2 donor0 246041 246044
#     PD-set1-2 donor1 246048 243003 
#     PD-set1-2 donor2 250118 250149
#     ...
# The input h5ad files look like:
#                              batch rep         set SubID_vs_w_gt SubID_vs_wo_gt
# AAACCCAAGCAATAAC  PD_set1_1_cDNA  1  PD_set1        246048         donor1
# AAACCCAAGCCTGAGA  PD_set1_1_cDNA  1  PD_set1      Negative       Negative
# AAACCCAAGCGCACAA  PD_set1_1_cDNA  1  PD_set1        donor2         donor0
# AAACCCAAGCGTTAGG  PD_set1_1_cDNA  1  PD_set1        246041         donor2
# AAACCCAAGGACATCG  PD_set1_1_cDNA  1  PD_set1      246041         donor2
# ...

# saved as /sc/arion/projects/CommonMind/pnm/AMP_PD/test_combo.h5ad, which looks like
#                              batch	rep	set	SubID_vs_w_gt	SubID_vs_wo_gt	donor_w_gt
# 	brain_reg_w_gt	donor_wo_gt	brain_reg_wo_gt
# AAACCCAAGCAATAAC	PD_set1_1_cDNA	1	PD_set1	246048	donor1	243003	BM-4	243003	BM-4
# AAACCCAAGCCTGAGA	PD_set1_1_cDNA	1	PD_set1	Negative	Negative				
# AAACCCAAGCGCACAA	PD_set1_1_cDNA	1	PD_set1	donor2	donor0				
# AAACCCAAGCGTTAGG	PD_set1_1_cDNA	1	PD_set1	246041	donor2	246044	BM-17	246044	BM-17
# AAACCCAAGGACATCG	PD_set1_1_cDNA	1	PD_set1	246041	donor2	246044	BM-17	246044	BM-17



# python3 create_aggregate_matrix.py /sc/arion/projects/CommonMind/pnm/AMP_PD/final_count_matrix/test/ \
#     /sc/arion/projects/psychAD/pnm/Hs_allchr_MT.txt \
#     /sc/arion/projects/CommonMind/pnm/AMP_PD/test_combo.h5ad \
#     --conversion_file conversion_vir_w_gt2.tsv \
#     --conversion_file_cols Sample,Donor_name,donor,Final_breg \
#     --suffix w_gt --h5ad_col ',brain_reg_' \
#     --conversion_file conversion_vir_first_run2.tsv \
#     --conversion_file_cols Sample,Donor_name,donor,Final_breg \
#     --suffix wo_gt --h5ad_col ',brain_reg_' \
#     --column_prefix SubID_vs_ -g GRCh38


import pegasus as pg, pegasusio as io, pandas as pd
import os, glob2, logging, re, argparse
from time import sleep
# import yaml


def get_names_anno(gid: str, col_n:str, df: pd.DataFrame):
    val=df.loc[df["gene_id"] == gid, col_n].values[0]
    if isinstance(val, str):
        return val
    elif isinstance(val, int):
        return str(val)
    else:
        return ""


def get_info(fn, ser1, conv_df):
    don_l = []
    breg_l = []
    for j in ser1:
        if j == "Negative" or j == "Doublet":
            don_l.append("")
            breg_l.append("")
            continue
        # if f == 'wo_gt':
        try:
            don_l.append(conv_df.loc[(conv_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
                & (conv_df["Donor_name"] == j), "donor"].values[0])
        except:
            don_l.append("")

        try:
            breg_l.append(conv_df.loc[(conv_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
                & (conv_df["Donor_name"] == j), "Final_breg"].values[0])
        except:
            breg_l.append("")
        # else:
        #     don_l.append(w_gt_df.loc[(w_gt_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
        #         & (w_gt_df["Donor_name"] == j), "donor"].values[0])
        #     breg_l.append(w_gt_df.loc[(w_gt_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
        #         & (w_gt_df["Donor_name"] == j), "Final_breg"].values[0])
            
    return don_l, breg_l


# def get_new_name(ser1):
#     ret_vals=[]
#     for x in ser1:
#         g_name = data.var.loc["gene_name"][data.var_names == x].values[0]
#         if (g_name).startswith('MT-'):
#                 ret_vals.append('MT-' + x.replace('_index', ''))
#         elif g_name != "":
#                 ret_vals.append(g_name+ '_' + re.sub(r"ENSG([0-9]+)_index", r"\1",  x))
#         else:
#                 ret_vals.append(x.replace('_index', ''))
# #         if g_name == "":
# #             ret_vals.append(x.replace('_index', ''))
#     return ret_vals


def get_new_name(x):
    if x[1].startswith('MT-'):
        return 'MT-' + x[0].replace('_index', '')
    else:
        return x[0].replace('_index', '')
    

def process_h5ad(fname, df_dict, inp_pref, 
    cnv_cols, h5ad_cols,
    gen="GRCh38"):
    data = pg.read_input(fname, genome=gen)
    chann = '-'.join(os.path.basename(fname).split('-')[:-1])
    for i, j in df_dict.items():
        data.obs['donor_' + i], data.obs['brain_reg_' + i] = get_info(chann, data.obs[inp_pref + i], j)
        data.obs['brain_reg_' + i] = data.obs['brain_reg_' + i].apply(lambda x: "GPi" if x == 'GPI' else x)

    return data


def get_argument_parser():
        """Generate and return argument parser."""

        #Parse Command-Line arguments
        parser = argparse.ArgumentParser(description="Demultiplex pools "
        "(supports hashsolo and vireo)"
        )
        parser.add_argument('count_mat_dir', help="Path to directory "
        "containing all h5ad files."
        )
        parser.add_argument('gene_file', help="Path to file that "
        "contains gene names and ids for annotation (tab-separated txt file)"
        )
        parser.add_argument('output_file', help="Path to store the final "
        "count matrix(h5ad)"
        )
        parser.add_argument('-c', '--conversion_file', action='append', 
        help="List of conversion files. Example [file1 [file2...]]"
        )
        parser.add_argument('--pool_col', action='append', help="Column "
        "(per conversion file), in the conversion files, containing "
        "pool names file"
        )
        parser.add_argument('--from_don_col', action='append', 
        help="Column (per conversion file), in the conversion files, "
        "containing donor names that need conversion"
        ""
        )
        parser.add_argument('--to_don_col', action='append', 
        help="Column (per conversion file), in the conversion files, "
        "containing the 'new' donor names"
        ""
        )
        parser.add_argument('--extra_anno_cols', 
        action='append', help="List of columns (per conversion file), "
        "in the conversion files, containing extra annotations. Example"
        " [annoA1,annoA2.. [annoB1,annoB2.. ]]", default=None, 
        nargs='?', const=None,
        )
        parser.add_argument('-s', '--suffix', required=True,
        action='append', help="List of suffixes for donor column names"
        " in the output file corresponding to the conversion files. "
        "Similar suffixes are expected for input file columns. Example"
        " [suffix1 [suffix2..]]"
        )
        parser.add_argument('--h5ad_extra_anno',
        action='append', help="List of column names in output file "
        "corresponding extra annotations (commma separated) present "
        "in conversion files [name1,name2.. [nameA,nameB]]. "
        "Should be provided for each extra annotation column. When "
        "no change is required use ''.", default=None, nargs='?', 
        const=None,
        )
        parser.add_argument('--column_prefix', 
        help="Column prefix used in input h5ads that contain"
        "IDs, which need to be converted. When multiple columns are "
        "needed this value is assumed to be the same for those columns"
        )
        parser.add_argument('-g', '--genome', help="Genome. Can be "
        "GRCh38, etc.", default='GRCh38',
        )

        return parser


def main():
    """Main entry point"""

    parser = get_argument_parser()
    args = parser.parse_args()

#     with open("../new_config.yaml") as fout:
#     sample_dict = yaml.load(fout, Loader=yaml.SafeLoader)

    op_cols = args.suffix
    inp_col_pref = args.column_prefix
    # inp_ext_anno = {} (
    #     args.extra_anno_cols.split(',')
    #     if args.extra_anno_cols is not None
    #     else args.extra_anno_cols
    #     )
    
    cnv_cols = dict.fromkeys(op_cols)
    per_file_key_list = ['file', 'inp_cols', 'conv_cols']
    # h5ad_cols = {} (
    #     args.h5ad_extra_anno.split(',')
    #     if args.h5ad_extra_anno is not None
    #     else args.h5ad_extra_anno
    #     )
    
    # Sanity checks
    assert len(args.conversion_file) == len(args.suffix), \
    "Number of column names in final h5ad and number of conversion " \
    "files are not same"
    assert len(args.extra_anno_cols) == len(args.conversion_file), \
    "Number of extra annotation columns provided are not same as the " \
    "number of conversion files provided"
    assert len(args.h5ad_extra_anno) == len(args.conversion_file), \
    "Number of extra annotation columns, for the final h5ad file, " \
    " provided are not same as the number of conversion files provided"
    # assert len(inp_ext_anno) == len(h5ad_cols), \
    # "Number of column names in final h5ad and number of conversion " \
    # "files are not same"

    mat_dir = os.path.join(args.count_mat_dir, "**/*.h5ad")
    all_f = glob2.glob(mat_dir)
    all_f.sort()
    all_f = [f for f in all_f if os.path.basename(f).split('_')[0] not in \
        ['PD-Set45-E1-HTO', 'PD-Set37-E2-HTO', 'PD-Set21-C2-HTO', 'PD-Set81-E1-HTO']]
    
    print(len(all_f))

    t2g = pd.read_csv(args.gene_file, skiprows=1, 
        names=["gene_id", "gene_name", "gene_start", "gene_end", "chr", "gene_type"], 
        sep="\t")
    t2g.index = t2g.gene_id

    # conv_df={}
    # for i, j in enumerate(args.conversion_file):
    #     temp_df = pd.read_csv(j, sep="\t")
    #     temp_df.fillna('', inplace=True)
    #     conv_df[op_cols[i]] = temp_df
    #     cnv_cols[op_cols[i]] = args.conversion_file_cols[i].split(',')

    # Create dict with suffixes for each conversion file
    conv_df = dict.fromkeys(op_cols)
    for i, j in conv_df.items():
        conv_df[i] = dict.fromkeys(per_file_key_list)

    for i, j in enumerate(args.conversion_file):
        temp_df = pd.read_csv(j, sep="\t")
        temp_df.fillna('', inplace=True)
        # conv_df[op_cols[i]] = temp_df
        conv_df[op_cols[i]] = dict.fromkeys(op_cols)
        conv_df[op_cols[i]]['file'] = j
        conv_df[op_cols[i]]['inp_cols'] = (
            args.extra_anno_cols[i].split(',')
            if args.extra_anno_cols is not None \
            and args.extra_anno_cols[i] is not None
            else None
            )
        conv_df[op_cols[i]]['conv_cols'] = (
            args.h5ad_extra_anno[i].split(',')
            if args.h5ad_extra_anno is not None \
            and args.h5ad_extra_anno[i] is not None
            else None
            )

        # cnv_cols[op_cols[i]] = args.conversion_file_cols[i].split(',')

    file_dict = dict.fromkeys(['Sample', 'Object'])
    file_dict['Sample'] = [ '_'.join(os.path.basename(f).split('-')[1:-1]) for f in all_f]
    file_dict['Object'] = [ 
        process_h5ad(f, conv_df, inp_col_pref, conv_df, args.genome) \
        for f in all_f
    ]
    # file_dict['Object'] = [ 
    #     process_h5ad(f, conv_df, inp_col_pref, args.genome) \
    #     for f in all_f
    # ]

    data = pg.aggregate_matrices(file_dict, default_ref=args.genome)
    data.var['gene_name'] = [ get_names_anno(g.split('_')[0], 'gene_name', t2g) for g in data.var.index.tolist()]
    data.var['chr'] = [ get_names_anno(g.split('_')[0], 'chr', t2g) for g in data.var.index.tolist()]
    data.var['gene_type'] = [ get_names_anno(g.split('_')[0], 'gene_type', t2g) for g in data.var.index.tolist()]

    data.var_names = data.var["gene_name"].reset_index().apply(get_new_name, axis=1).values

    io.write_output(data, args.output_file)


if __name__ == '__main__':

    main()
    sleep(20)