#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import glob2, os, re, argparse
from time import sleep

 


def read_files_ext(fname) -> pd.DataFrame :
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 
    if fname.endswith('.csv'):
        return pd.read_csv(fname)
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t')
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")



def set_don_ids(x) -> str:
    if x == 'doublet':
        return 'Doublet'
    elif x == 'unassigned':
        return 'Negative'
    elif x.startswith('donor'):
        return x
    else:
        return '_'.join(x.split('_')[1:])


def get_don_ids(x, t_df) -> str:
    try:
        return t_df.loc[t_df["primary_genotype"] == x, "SubID"].values[0]
    except:
        return x


def ret_subj_ids(ser, t_df) -> pd.DataFrame:
    headers = ["Subj_ID", "prob_max", "prob_doublet"]
    ret_df_l = []
    for x in ser:
        try:
            ret_df_l.append(t_df.loc[t_df["cell"].str.strip() == x.strip(), headers].values.flatten().tolist())
        except:
            ret_df_l.append(["Not Present", "NA", "NA"])

    return pd.DataFrame(ret_df_l, columns=headers)




sc.settings.set_figure_params(dpi_save=400, format='png', color_map = 'viridis_r')
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()


def get_argument_parser():
    """Generate and return argument parser."""
        
    #Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Demultiplex sample based on vireoSNP produced output")

    parser.add_argument('h5ad_file', help="h5ad file")
    parser.add_argument('donors_file', help="donor_ids.tsv")
    parser.add_argument('converter_file', help="File need for converting genotype IDs to Subject IDs")
    parser.add_argument('output_file', help="output_h5ad file")

    return parser


def main():

    parser = get_argument_parser()
    args = parser.parse_args()


    # Storing parsed inputs
    ann = ad.read(args.h5ad_file)
    vir_class = read_files_ext(args.donors_file)
    conv_df = pd.read_csv(args.converter_file)
    op = args.output_file

    # Create necessary folders
    if not os.path.isdir(op.replace('/' + os.path.basename(op), '')):
        os.makedirs(op.replace(os.path.basename(op), ''))


    # vir_class.rename(columns={"cell":"barcodes", "donor_id":"Subj_ID"}, inplace=True, errors="raise")
    vir_class["donor_id"] = vir_class["donor_id"].apply(set_don_ids)
    conv_df = conv_df.loc[conv_df['primary_genotype'].isin(vir_class['donor_id'].unique()), ["SubID", "primary_genotype"]]
    vir_class['Subj_ID'] = vir_class['donor_id'].apply(get_don_ids, args=(conv_df,))
    del vir_class['donor_id']

    # get_df = adata.obs_names.to_series().apply(ret_subj_ids, args=(vir_class,))
    get_df = ret_subj_ids(ann.obs_names.to_list(), vir_class)
    ann.obs['SubID_vS'] = get_df["Subj_ID"].to_list()
    ann.obs['max_prob'] = get_df["prob_max"].to_list()
    ann.obs['doublet_prob'] = get_df["prob_doublet"].to_list()

    ann.write(op)


if __name__ == '__main__':
    main()
    sleep(60)
