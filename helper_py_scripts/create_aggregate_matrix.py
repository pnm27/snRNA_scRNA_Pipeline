#!/usr/bin/env python3


import pegasus as pg, pegasusio as io, pandas as pd, os, glob2, logging, yaml, argparse


def get_names_anno(gid, col):
    val=t2g.loc[t2g["gene_id"] == gid, col].values[0]
    if isinstance(val, str):
            return val
    elif isinstance(val, int):
            return str(val)
    else:
            return ""


def get_info(fn, ser1, f=None):
    don_l = []
    breg_l = []
    for j in ser1:
            if j == "Negative" or j == "Doublet":
                    don_l.append("")
                    breg_l.append("")
                    continue
            if f == 'wo_gt':
                    don_l.append(wo_gt_df.loc[(wo_gt_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
                    & (wo_gt_df["Donor_name"] == j), "donor"].values[0])
                    breg_l.append(wo_gt_df.loc[(wo_gt_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
                    & (wo_gt_df["Donor_name"] == j), "Final_breg"].values[0])
            else:
                    don_l.append(w_gt_df.loc[(w_gt_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
                    & (w_gt_df["Donor_name"] == j), "donor"].values[0])
                    breg_l.append(w_gt_df.loc[(w_gt_df["Sample"].str.lower().str.contains(fn.lower(), na=False)) \
                    & (w_gt_df["Donor_name"] == j), "Final_breg"].values[0])
    return don_l, breg_l


def process_h5ad(fname, gen="GRCh38"):
    data = pg.read_input(fname, genome=gen)
    chann = '-'.join(os.path.basename(fname).split('-')[:-1])
    data.obs['donor_wo_gt'], data.obs['brain_reg_wo_gt'] = get_info(chann, data.obs['SubID_vs_wo_gt'], "wo_gt")
    data.obs['donor_w_gt'], data.obs['brain_reg_w_gt'] = get_info(chann, data.obs['SubID_vs_w_gt'])
    data.obs['brain_reg_w_gt'] = data.obs['brain_reg_w_gt'].apply(lambda x: "GPi" if x == 'GPI' else x)
    data.obs['brain_reg_wo_gt'] = data.obs['brain_reg_wo_gt'].apply(lambda x: "GPi" if x == 'GPI' else x)
    return data


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
    "Default: None", default=None,
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


    with open("../new_config.yaml") as fout:
        sample_dict = yaml.load(fout, Loader=yaml.SafeLoader)


    mat_dir = "/sc/arion/projects/CommonMind/pnm/AMP_PD/final_count_matrix/vireoSNP_1kGP_isec/**/*.h5ad"
    all_f = glob2.glob(mat_dir)
    all_f.sort()
    print(len(all_f))
    all_f = [f for f in all_f if os.path.basename(f).split('_')[0] not in ['PD-Set45-E1-HTO', 'PD-Set37-E2-HTO', 'PD-Set21-C2-HTO', 'PD-Set81-E1-HTO']]
    print(len(all_f))

    gene_info_dir = "/sc/arion/projects/psychAD/pnm/Hs_allchr_MT-DL.txt"
    t2g = pd.read_csv(gene_info_dir, skiprows=1, names=["gene_id", "gene_name", "gene_start", "gene_end", "chr", "gene_type"], sep="\t")
    t2g.index = t2g.gene_id

    wo_gt_df = pd.read_csv("conversion_vir_first_run2.tsv", sep="\t")
    w_gt_df = pd.read_csv("conversion_vir_w_gt2.tsv", sep="\t")
    wo_gt_df.fillna('', inplace=True)
    w_gt_df.fillna('', inplace=True)

    file_dict = dict.fromkeys(['Sample', 'Object'])
    file_dict['Sample'] = [ '_'.join(os.path.basename(f).split('-')[1:-1]) for f in all_f]
    file_dict['Object'] = [ process_h5ad(f) for f in all_f]

    data = pg.aggregate_matrices(file_dict, default_ref="GRCh38")

    io.write_output(data, 'all_cells_fixed_annoSwap_v2.h5ad')