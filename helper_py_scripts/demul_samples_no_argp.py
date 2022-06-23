#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

# Solo didn't run through scvi, scvi-tools nor scanpy.external
# Only this seems to work
from typing import Union # Need verion > 3.5
import anndata as ad
import scanpy as sc, pandas as pd, numpy as np
import glob2, os, re, sys
from collections import Counter
from openpyxl import load_workbook
from collections import defaultdict, OrderedDict as ord_dict
import datetime
import functools
from time import sleep
 

assert sys.version_info >= (3, 5), "This script needs python version >= 3.5!"


# Basic helper functions-------------------------------------------------------------------------------------------
def read_files_ext(fname) -> pd.DataFrame :
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 
    if fname.endswith('.csv'):
        return pd.read_csv(fname)
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t')
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")


# Simplify these 'parse' functions
def parse_HTO(wet_lab_df, col_val, fname, s_name, hs=None) -> list:
    sub = wet_lab_df[col_val]
    test_len = len(sub)
    # No command-line params and inference that all HTOs are present in one row separated by ","
    if hs == None and test_len == 1 and sub.str.count(',').values[0] > 1:
        return sub.values[0].split(',')
    # No command-line params and inference that all HTOs are present in one row separated by whitespaces    
    elif hs == None and test_len == 1 and len(sub.split()) > 1:
        return sub.values[0].split(',')
    elif hs == None and test_len > 1:
        return sub.tolist()
    elif hs != None and test_len == 1 and sub.str.count(hs).values[0] > 1:
        return sub.values[0].split(hs)
    elif hs != None and test_len > 1:
        raise ValueError(f"After subsetting sample {s_name} from the wet lab file {fname}, there are multiple rows ({test_len}) of HTOs for this sample while a separator value is also provided.")
    elif hs != None and test_len == 1 and sub.str.count(hs).values[0] == 1:
        raise ValueError(f"Either the given separator {hs} is wrong or the sample {s_name} has incomplete HTO values in the wet lab file {fname}")
    else:
        raise ValueError(f"Something is wrong with the given input(s):\n\twet lab file: {fname}\n\tsample: {s_name}\n\tHTO-separator: {hs}")


# Assumes similar construct like the HTOs
def parse_subids(wet_lab_df, col_val, fname, s_name, hs=None) -> list:
    sub = wet_lab_df[col_val]
    test_len = len(sub)
    # No command-line params and inference that all HTOs are present in one row separated by ","
    if hs == None and test_len == 1 and sub.str.count(',').values[0] > 1:
        return sub.values[0].split(',')
    # No command-line params and inference that all HTOs are present in one row separated by whitespaces    
    elif hs == None and test_len == 1 and len(sub.split()) > 1:
        return sub.values[0].split(',')
    elif hs == None and test_len > 1:
        return sub.tolist()
    elif hs != None and test_len == 1 and sub.str.count(hs).values[0] > 1:
        return sub.values[0].split(hs)
    elif hs != None and test_len > 1:
        raise ValueError(f"After subsetting sample {s_name} from the wet lab file {fname}, there are multiple rows ({test_len}) of HTOs for this sample while a separator value is also provided.")
    elif hs != None and test_len == 1 and sub.str.count(hs).values[0] == 1:
        raise ValueError(f"Either the given separator {hs} is wrong or the sample {s_name} has incomplete donor ids in the wet lab file {fname}")
    else:
        raise ValueError(f"Something is wrong with the given input(s):\n\twet lab file: {fname}\n\tsample: {s_name}\n\tHTO-separator: {hs}")


def redo_test(inp_key, rem_last_n_chars=0) -> Union[str, None]:
    if inp_key in snakemake.input.keys() and os.path.isfile(snakemake.input[inp_key]) and rem_last_n_chars==0:
        return snakemake.input[inp_key]
    elif inp_key in snakemake.input.keys() and os.path.isfile(snakemake.input[inp_key]) and rem_last_n_chars>0:
        return snakemake.input[inp_key][:-rem_last_n_chars]
    elif inp_key in snakemake.input.keys() and not os.path.isfile(snakemake.input[inp_key]):
        raise ValueError(f"The file {snakemake.input[inp_key]} does not exist!")
    else:
        return None


# For calico_solo------------------------------------------------------------------------
def ret_htos_calico_solo(bcs, df_s):

            # List of htos from the wet lab spreadsheet
            hto_l = parse_HTO(df_s, cols[1], snakemake.params.samples_info, samp, snakemake.params.hto_sep)
            # List of subIDs from the wet lab spreadsheet
            subid_l = parse_subids(df_s, cols[1], snakemake.params.samples_info, samp, snakemake.params.hto_sep)

            # List of barcodes
            barc_l = []
            # SubID from Shan's csv file
            ret_samp = []
            # List of HTOs as HTO1, HTO2, etc
            hash_n = []
            # Doublets' count
            doublet_n = 0
            # Negatives' count
            negative_n = 0

            for bc in bcs:
                if bc in dem_cs.obs_names:
                    hto_n = dem_cs.obs[dem_cs.obs_names == bc].Classification.values[0]
                    
                   #samp_n.append(df_shan[(df_shan["cDNA_ID"] == b) & (df_shan["hashtag" == hto_n]), "sample_ID"].values[0])
                    if hto_n == 'Doublet':
                        barc_l.append(bc)
                        hash_n.append(hto_n)
                        ret_samp.append(hto_n)
                        doublet_n += 1

                    elif hto_n == 'Negative':
                        barc_l.append(bc)
                        hash_n.append(hto_n)
                        ret_samp.append(hto_n)
                        negative_n += 1

                    else:
                        barc_l.append(bc)
                        hash_n.append(hto_n)
                        try:
                            ret_samp.append(subid_l[hto_l.index(hto_n)])
                        except:
                            ret_samp.append(hto_n)


            
            ser_s = pd.DataFrame({'Sample':ret_samp, 'HTO':hash_n}, index=barc_l)

            return [ser_s, doublet_n, negative_n] #zip(barc_l, samp_n)



        def ret_samp_names(y, df_info):
            #print(len(y), len(bc))
            if y in df_info.index:
                return df_info.loc[df_info.index == y, 'Sample'].values[0]
            
            else:
                return "Not Present"



        def ret_hto_number(y, df_info):
            #print(len(y), len(bc))
            if y in df_info.index:
                return df_info.loc[df_info.index == y, 'HTO'].values[0]

            else:
                return "Not Present"


# For adding vireoSNP classfication
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


#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Demultiplex sample based on hahsolo produced output")

parser.add_argument('matrix_file', help="Path to matrix.mtx.gz")
parser.add_argument('count_matrix', help="Path to store the final count matrix(h5ad)")
parser.add_argument('demux_info', help="Path to store demultiplexing info (tab-separated txt file)")
parser.add_argument('gene_info_file', help="Path to file that contains gene names and ids for annotation (tab-separated txt file)")
parser.add_argument('wet_lab_file', help="Path to file that contains HTO info for each set (either csv or tsv file)")


# Optional parameters
parser.add_argument('--calico_solo', dest='hashsolo_out', help="Path to cached output of hashsolo(h5ad)")
parser.add_argument('--vireo_out', help="Path to cached output of hashsolo(h5ad)")
parser.add_argument('-m', '--max_mito', type=int, help="Max mitochondrial genes(in percent) per cell. Default: 5", default=5)
parser.add_argument('-g', '--min_genes', type=int, help="Min #genes per cell. Default: 1000", default=1000)
parser.add_argument('-c', '--min_cells', type=int, help="Min #cells expressing a gene for it to pass the filter. Default: 10", default=10)
parser.add_argument('--columns', nargs=4, help="List of column names RESPECTIVELY to sample_ID (as present in the spreadsheets - identifies pooled samples), \
    HTO numbers, it's associated barcodes, and Donors/SubIDs (contains each multiplexed donors).", metavar=('sample_ID', 'HTO_name', 'HTO_barcode', 'Sub_ID'),
     default=['unique_sample_ID', 'hashtag', 'ab_barcode', 'SubID'])
parser.add_argument('-s', '--sample_name', help="Name of the sample. Check whether this name coincides with that in this script as well as the one in the wet_lab_file")
parser.add_argument('--hto_sep', help="If, per each sample in the wet lab file (6th positional argument to this script), HTOs are all present in one row separated by some SEP then specify it here. Default: ' '", default=' ')


args = parser.parse_args()

# Adding demultiplex info or creating "new" final count matrix 
redo = redo_test('prev_count_mtx')

# If 'redoing' demultiplexing then whether vireo output is to be added or calico_solo output is to be added to an existing final count matrix file
add_calico = redo_test('hashsolo_out') if redo is not None else None
add_vireo = redo_test('vireo_out') if redo is not None else None


# If creating "new" final count matrix then selecting both or one
starsolo_mat = redo_test('starsolo_out', 13) if redo is None else None #args.matrix_file[:-13]
calico_demux = redo_test('hashsolo_out') if redo is None else None
vireo_demux = redo_test('vireo_out') if redo is None else None

# Can't have both 'starsolo_out' and 'prev_count_mtx' as None or not None
assert redo is not None and starsolo_mat is not None, "Previously run count matrix file is provided along with STARsolo's output. Please give one depending on whether redoing or first demultiplex run!"
assert redo is None and starsolo_mat is None, "Neither previously run count matrix file is provided nor the STARsolo's output. Please give ONE depending on whether redoing or first demultiplex run!"

# Store output_file and Create necessary folders
op = snakemake.output[0] #args.count_matrix
if not os.path.isdir(op.replace('/' + os.path.basename(op), '')):
    os.makedirs(op.replace(os.path.basename(op), ''))

# Initial run of demultiplexing
if redo is None:
    # Both calico_solo and vireo can't be None
    if calico_demux is None and vireo_demux is None:
        raise ValueError("For the inital demultiplexing run for STARsolo's output provide either a VALID calico_solo h5ad or vireoSNP's donor file")

    cols = snakemake.params.cols #args.columns

    # Batch info 
    batch=snakemake.params.sample_name.replace('-', '_')+'_cDNA'

    # Prepare for Extra Information
    r_num=int(os.path.basename(op).split('_')[0].replace('round', ''))
    preparer_rep=batch.split('_')[2]
    replicate=preparer_rep[1]
    preparer=preparer_rep[0]

    # batch for wet lab file 
    samp=batch

    # Parameters for filtering
    max_mito=snakemake.params.mito
    min_genes=snakemake.params.min_genes
    min_cells=snakemake.params.min_cells

    t2g = pd.read_csv(snakemake.params.genes_info, skiprows=1, usecols=range(2),names=["gene_id", "gene_name"], sep="\t")
    t2g.index = t2g.gene_id
    t2g = t2g.loc[~t2g.index.duplicated(keep='first')]


    # Wet Lab file, Filter wet lab file's columns, if needed
    df = read_files_ext(snakemake.params.samples_info)
    if df.loc[df[cols[0]] == samp].empty:
        raise ValueError(f"Check dtypes!\nSample (variable name 'var', data type {type(samp)}, with value {samp} ) couldn't be subset from the wet lab file.\nData types for the wet lab file:\n{df.dtypes}")
    else:
        df = df.loc[df[cols[0]] == samp]


    # Process STARsolo output----------------------------------------------------------------------------------------------------------------------------
    ct = datetime.datetime.now()
    print(f"Processing STARsolo's output at: {ct}")
    try:
        adata = sc.read_10x_mtx(starsolo_mat, make_unique=True, var_names= "gene_ids", cache=True)
    except:
        e = sys.exc_info()[0]
        print(f"Error encountered while loading the mtx files!\nError message: {e}")

    solo_run_info = []
    print(adata)
    solo_run_info.append(( 'Started with cells', adata.n_obs))
    solo_run_info.append(( 'Started with genes', adata.n_vars))
    adata.var_names_make_unique()
    adata.var["gene_id"] = adata.var.index.values
    adata.var["gene_name"] = adata.var.gene_id.map(t2g["gene_name"])
    adata.var_names = adata.var_names.to_series().map(lambda x: x + '_index')
    adata = adata[:, pd.notna(adata.var["gene_name"])]
    adata.var_names_make_unique()
    adata.X = adata.X.astype('float64')
    solo_run_info.append(( 'Unique genes', adata.n_vars))
    print(adata)

    # Filter data using cell level metrics
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    solo_run_info.append(( 'min #genes expressed per cell', min_genes))
    solo_run_info.append(( 'min #cells expressing per gene', min_cells))
    solo_run_info.append(( 'Retained cells after previous filter', adata.n_obs))
    solo_run_info.append(( 'Retained genes after previous filter', adata.n_vars))

    # Filter data wrt mito content
    adata.var["mito"] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
    adata = adata[adata.obs["pct_counts_mito"]< max_mito, :]
    solo_run_info.append(( 'max percent mito content per cell', max_mito))
    solo_run_info.append(( 'cells with low mito percent', adata.n_obs))
    print(adata)

    # Conditional execution dependent on demultiplexing methods selected
    if calico_demux is not None and vireo_demux is not None:
        pass # to do

    elif calico_demux is not None:

        # Load hashsolo/calico_solo output (h5ad)
        try:
            dem_cs = ad.read(calico_demux)
        except:
            e = sys.exc_info()[0]
            print(f"Error encountered while loading the h5ad file!\nError message: {e}")

        print("Successfully loaded calico_solo output!")

        # Demultiplex and assign samples----------------------------------------------------------

        # hto_tags_cs, hto_tags_ms, and hto_tags_hd contain (in sequence): pd DF with barcodes as index. SubID and HTO number (HTO1, HTO2, etc) as columns
        #, no. of doublet cells, no. of negative cells]
        ct = datetime.datetime.now()
        print(f"Starting: Calculating demultiplexing info for hashsolo/calico solo at: {ct}")
        hto_tags_cs = ret_htos_calico_solo(adata.obs_names.to_list(), df)
        ct = datetime.datetime.now()
        print(f"Finished: Calculating demultiplexing info for hashsolo/calico solo at: {ct}")


        # add few more annotations
        adata.obs['batch'] = batch
        adata.obs['round_num'] = r_num
        adata.obs['prep'] = preparer
        adata.obs['rep'] = replicate
        adata.obs['set'] = batch[:-6]

        # Create obs columns in adata to represent the SubID as assigned by calico solo and its associated hastag number
        SubID_cs = adata.obs_names.to_series().apply(ret_samp_names, args=(hto_tags_cs[0], ))
        hasht_n_cs = adata.obs_names.to_series().apply(ret_hto_number, args=(hto_tags_cs[0], ))
        adata.obs['SubID_cs'] = SubID_cs
        adata.obs['HTO_n_cs'] = hasht_n_cs

        # Save doublets and negatives info from calico solo
        solo_run_info.append(('Doublets #cells_cs', hto_tags_cs[1]))
        solo_run_info.append(('Negative #cells_cs', hto_tags_cs[2]))


        # Remaining cells after demultiplexing (for each demux method)
        rem_cells_cs = adata.obs.SubID_cs.value_counts()
        prop_dict = rem_cells_cs[(rem_cells_cs.index != "Doublet") & (rem_cells_cs.index != "Negative") & (rem_cells_cs.index != "Not Present")]
        od = ord_dict(sorted(prop_dict.items()))
        a = ""
        for k, v in od.items():
            b = str(k) + ": " + str(v) + ", "
            a += b

        solo_run_info.append(( 'After demultiplexing #cells_cs', a[:-1] ))


        ct = datetime.datetime.now()
        print(f"Saving All info as a tsv file and also the h5ad files: {ct}")
        solo_run_df = pd.DataFrame(solo_run_info, columns=['Observations', 'Vals'])
        solo_run_df.to_csv(snakemake.output[1], sep = "\t", index=False)

        # If Subject IDs aren't 'string' then convert them
        adata.obs['SubID_cs']=adata.obs['SubID_cs'].apply(str)

        adata.write(op)

        ct = datetime.datetime.now()
        print(f"Finished: Processing Sample {batch} at: {ct}")

    # To do
    elif vireo_demux is not None:
        print("Starting demultiplexing through vireoSNP's output")
        # Demultiplex and assign samples----------------------------------------------------------

        # hto_tags_cs, hto_tags_ms, and hto_tags_hd contain (in sequence): pd DF with barcodes as index. SubID and HTO number (HTO1, HTO2, etc) as columns
        #, no. of doublet cells, no. of negative cells]
        ct = datetime.datetime.now()
        print(f"Starting: Calculating demultiplexing info for hashsolo/calico solo at: {ct}")
        hto_tags_cs = ret_htos_calico_solo(adata.obs_names.to_list(), df)
        ct = datetime.datetime.now()
        print(f"Finished: Calculating demultiplexing info for hashsolo/calico solo at: {ct}")


        # add few more annotations
        adata.obs['batch'] = batch
        adata.obs['round_num'] = r_num
        adata.obs['prep'] = preparer
        adata.obs['rep'] = replicate
        adata.obs['set'] = batch[:-6]

        # Create obs columns in adata to represent the SubID as assigned by calico solo and its associated hastag number
        SubID_cs = adata.obs_names.to_series().apply(ret_samp_names, args=(hto_tags_cs[0], ))
        hasht_n_cs = adata.obs_names.to_series().apply(ret_hto_number, args=(hto_tags_cs[0], ))
        adata.obs['SubID_cs'] = SubID_cs
        adata.obs['HTO_n_cs'] = hasht_n_cs

        # Save doublets and negatives info from calico solo
        solo_run_info.append(('Doublets #cells_cs', hto_tags_cs[1]))
        solo_run_info.append(('Negative #cells_cs', hto_tags_cs[2]))


        # Remaining cells after demultiplexing (for each demux method)
        rem_cells_cs = adata.obs.SubID_cs.value_counts()
        prop_dict = rem_cells_cs[(rem_cells_cs.index != "Doublet") & (rem_cells_cs.index != "Negative") & (rem_cells_cs.index != "Not Present")]
        od = ord_dict(sorted(prop_dict.items()))
        a = ""
        for k, v in od.items():
            b = str(k) + ": " + str(v) + ", "
            a += b

        solo_run_info.append(( 'After demultiplexing #cells_cs', a[:-1] ))


        ct = datetime.datetime.now()
        print(f"Saving All info as a tsv file and also the h5ad files: {ct}")
        solo_run_df = pd.DataFrame(solo_run_info, columns=['Observations', 'Vals'])
        solo_run_df.to_csv(snakemake.output[1], sep = "\t", index=False)

        # If Subject IDs aren't 'string' then convert them
        adata.obs['SubID_cs']=adata.obs['SubID_cs'].apply(str)

        adata.write(op)

        ct = datetime.datetime.now()
        print(f"Finished: Processing Sample {batch} at: {ct}")
    else:
        raise ValueError("Unexpected condition encountered! (bug: 1)") # Should never be executed

elif redo is not None:
    try:
        ann = ad.read(snakemake.input.h5ad_file)
    except:
        e = sys.exc_info()[0]
        print(f"Error encountered while loading the h5ad file (previously completed demultiplex run)!\nError message: {e}")    

    if add_vireo is not None:
        # Storing parsed inputs        
        vir_class = read_files_ext(args.donors_file)
        conv_df = pd.read_csv(args.converter_file)

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
        sleep(60)

else:
    raise ValueError("Unexpected condition encountered! (bug: 2)") # Should never be executed



sleep(60)