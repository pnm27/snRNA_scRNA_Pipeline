#!/usr/bin/env python
# coding: utf-8


from __future__ import annotations
import pandas as pd, os, argparse, glob2, matplotlib, numpy as np
from adjustText import adjust_text # To avoid overlapping texts

# Anti-Grain Geometry, which is a non-GUI backend used for writing images to files, rather than displaying them on a screen
# BEFORE PYPLOT IS CALLED
matplotlib.use("Agg")

import matplotlib.pyplot as plt

# Multiprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def process_file(f, op_png_dir):

    # Mark up to 4 possible contenders to best possible (sample with highest %hets, followed by %homo):
    # a) top 2 samples with highest %hom
    # b) top 2 samples with highest %hets
    samp = os.path.basename(f).replace(".bamstat.txt", "").split('_')[0]
    don = '_'.join(os.path.basename(f).replace(".bamstat.txt", "").split('_')[1:])
    op_dir = os.path.join(op_png_dir, samp)
    os.makedirs(op_dir, exist_ok=True)

    op_n = os.path.join(op_dir, f"{don}.png")

    if os.path.isfile(op_n):
        return f"Skipping {f}"
    
    f1=pd.read_csv(f, sep=' ')
    # 1. Find the point closest to (1,1)
    f1["dist_to_11"] = np.sqrt((1 - f1["x"])**2 + (1 - f1["y"])**2)
    
    plot_samples = [] # To check if a sample is already plotted or not
    texts = [] # To align the texts in the plot
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    ax.scatter(f1['perc_het_consistent'], f1['perc_hom_consistent'])

    f1_sorted = f1.sort_values(
        ['perc_het_consistent', 'perc_hom_consistent'],
        ascending=False
    ).reset_index(drop=True)

    for i in [0, 1]:
        samp_id = f1_sorted['SampleID'][i]
        x, y = f1_sorted.iloc[i][['perc_het_consistent', 'perc_hom_consistent']]
        if not any(np.isnan((x, y))):
            plot_samples.append(samp_id)
            texts.append(ax.annotate(samp_id, (x, y)))

    f1_sorted = f1.sort_values(
        ['perc_hom_consistent', 'perc_het_consistent'],
        ascending=False
    ).reset_index(drop=True)

    for i in [0, 1]:
        samp_id = f1_sorted['SampleID'][i]
        x, y = f1_sorted.iloc[i][['perc_het_consistent', 'perc_hom_consistent']]
        if samp_id not in plot_samples and not any(np.isnan((x, y))):
            plot_samples.append(samp_id)
            texts.append(ax.annotate(samp_id, (x, y)))

    
    # if not an empty list of samples to plot
    # 1) adjust texts
    # 2) Add a table below the fig
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))
        plot_df = f1.loc[
            f1['SampleID'].isin(plot_samples), 
            ["SampleID", "n_geno_missing", 
             "perc_het_consistent", "perc_hom_consistent"]
        ].copy()
        plot_df = plot_df.rename(columns = {
            "n_geno_missing": "n_geno_miss", 
            "perc_het_consistent": "perc_het", 
            "perc_hom_consistent": "perc_hom"
        })
        table = ax.table(
            cellText=plot_df.to_numpy().tolist(), 
            loc='bottom', 
            bbox=[0, -0.2, 1, 0.1], 
            colLabels=plot_df.columns
        )
        
        # Set column header text to bold
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.get_text().set_weight('bold')
    

    fig.tight_layout()
    plt.savefig(op_n, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close(fig)
    return f"Finished plotting for: {f}"


def get_argument_parser():
    """Generate and return argument parser."""

    #Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Demultiplex pools "
    "(supports hashsolo and vireo). Note: poor cells are not removed but "
    "they aren't included while demultiplexing."
    )
    parser.add_argument('qtltools_dir', help="Path to matrix.mtx.gz or h5ad "
    "file. If an h5ad file is provided then it is expected that it "
    "has been already processed i.e. poor cells are already filtered out."
    )
    parser.add_argument('--qtltools_suffix', nargs='?', help="Column containing "
    "converted donor names in the converter file. As present in snakemake's "
    "config['identify_swaps']['mbv_suffix']. When flag is used but no value "
    "is provided = .bamstat.txt, no flag provided = .txt", 
    const=".bamstat.txt", default='.txt'
    )
    parser.add_argument('--vireo_out', help="Path to the directory "
    "containing vireo outputs."
    )
    parser.add_argument('--vireo_sample_prefix', help="If the sample "
    "names in vireo outputs are different than that of outputs in "
    "qtltools."
    )
    parser.add_argument('--threads', help="Number of threads/parallel "
    "workers to use.", type=int, default=4,
    )
    parser.add_argument('-s', '--sep', help="Separator used in qtltools "
    "filenames between donor and pool. DEFAULT = '_'. "
    "example: poolA_donorX.bamstat.txt, poolA_donorY.bamstat.txt "
    "then sep='_'. First part of the split is c", default='_',
    )
    return parser


def main():

    parser = get_argument_parser()
    args = parser.parse_args()

    qtltools_dir = args.qtltools_dir
    qtltools_suff = args.qtltools_suffix

    # all_files=glob2.glob("/sc/arion/projects/CommonMind/pnm/genesis_U24/A6_Ruzicka/qtltools_mbv/*txt")
    all_files = glob2.glob(os.path.join(qtltools_dir, f"*{qtltools_suff}"))
    # op_png_dir="/sc/arion/projects/CommonMind/pnm/genesis_U24/A6_Ruzicka/qtltools_mbv_pngs/"
    op_png_dir = (
        qtltools_dir[:-1] + '_pngs/'
        if qtltools_dir[-1] == '/' 
        else qtltools_dir + '_pngs/'
    )
    # vir_out_dir = "/sc/arion/projects/psychAD/pnm/A6_Ruzicka/demultiplex/vireoSNP/"
    vir_out_dir = args.vireo_out

    max_workers = args.threads # 14   # Adjust to your CPU cores

    os.makedirs(op_png_dir, exist_ok=True)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_file, f, op_png_dir)
            for f in sorted(all_files)
        ]

        for future in as_completed(futures):
            print(future.result())


    a=[]
    uniq_id = {} # Unique for vireo's Sample and seq type pair
    for f in sorted(all_files):
        f1=pd.read_csv(f, sep=' ')
        f1.sort_values(['perc_het_consistent', 'perc_hom_consistent'], ascending=False, inplace=True)
        f1.reset_index(inplace=True, drop=True)
        # samp = os.path.basename(f).replace(".bamstat.txt", "").split('_')[0]
        # don = '_'.join(os.path.basename(f).replace(".bamstat.txt", "").split('_')[1:])
        samp = os.path.basename(f).replace(qtltools_suff, "").split('_')[0]
        don = '_'.join(os.path.basename(f).replace(qtltools_suff, "").split('_')[1:])
        # seq_type = os.path.basename(f).replace(".bamstat.txt", "").split('_')[0]
        # don = '_'.join(os.path.basename(f).replace(".bamstat.txt", "").split('_')[1:])
        if len(uniq_id.keys()) == 0 or os.path.join(vir_out_dir, f"Sample_{samp}/summary.tsv") not in uniq_id:
            summ_f = os.path.join(vir_out_dir, f"Sample_{samp}/summary.tsv")
            temp_f = pd.read_csv(summ_f, sep='\t')
            uniq_id[summ_f] = temp_f
        else:
            temp_f= uniq_id[os.path.join(vir_out_dir, f"Sample_{samp}/summary.tsv")]
        try:
            cell_cts = temp_f.loc[temp_f["Var1"] == don, "Freq"].values[0]
        except:
            cell_cts = 0
        # status = f1['SampleID'][0] in pool_vcf_conv.loc[pool_vcf_conv['Pool'] == '-'.join(samp.split('-')[:-1]), 'Matched vcf ID'].values
        # status = "In pool" if status else "Not in pool"
        # a.append((samp, seq_type, don, f1['SampleID'][0], f1['perc_het_consistent'][0], f1['perc_hom_consistent'][0], cell_cts, status))
        a.append((samp, don, f1['SampleID'][0], f1['perc_het_consistent'][0], f1['perc_hom_consistent'][0], cell_cts))


    # a_df = pd.DataFrame(a, columns=["Sample", "Seq_type", "Donor_name", "donID", "perc_het_consistency", "perc_hom_consistency", "cell_counts", "Status"])
    a_df = pd.DataFrame(a, columns=["Sample", "Donor_name", "donID", "perc_het_consistency", "perc_hom_consistency", "cell_counts"])
    a_df.to_csv("/sc/arion/projects/CommonMind/pnm/genesis_U24/A6_Ruzicka/A6_Ruzicka_wo_gt_res.tsv", sep='\t', index=False)

    def ret_strings(x, sep=','):
        return sep.join(x)


    b=[]
    for f in sorted(all_files):
        f1=pd.read_csv(f, sep=' ')
        f1.sort_values(['perc_het_consistent', 'perc_hom_consistent'], ascending=False, inplace=True)
        f1.reset_index(inplace=True, drop=True)
        val1 = f1.iloc[0, 8] - f1.iloc[1, 8]
        val2 = f1['perc_het_consistent'].max() - f1['perc_het_consistent'].min()
        samp = os.path.basename(f).replace(".bamstat.txt", "").split('_')[0]
        don = '_'.join(os.path.basename(f).replace(".bamstat.txt", "").split('_')[1:])

        b.append((samp, don, f1['SampleID'][0], val1, val2))


    b_df = pd.DataFrame(b, columns=["Sample", "Donor_name", "donID", "diff_het_top2vals", "range_het_vals"])
    b_df.to_csv("/sc/arion/projects/CommonMind/pnm/genesis_U24/A6_Ruzicka/A6_Ruzicka_wo_gt_res2.tsv", sep='\t', index=False)

    all_files=glob2.glob(os.path.join(vir_out_dir, "**/summary.tsv"))
    c = []
    for f in sorted(all_files):
        samp = f.split('/')[-2].replace('Sample_', '')
        temp_f = pd.read_csv(f, sep='\t')
            
        t_dons_l = temp_f["Var1"].tolist()
        t_cc_l = temp_f["Freq"].tolist()
        dons=[]
        cells=[]
        doub_val=0
        neg_val=0
        for i, j in enumerate(t_dons_l):
            if j != 'doublet' and j != 'unassigned':
                dons.append(j)
                cells.append(str(t_cc_l[i]))
            elif j == 'doublet':
                doub_val = t_cc_l[i]
            elif j == 'unassigned':
                neg_val = t_cc_l[i]
                
        c.append((samp, ret_strings(dons), ret_strings(cells), neg_val, doub_val ))


    c_df = pd.DataFrame(c, columns=["Sample", "donors", "cell_counts", "n_negs", "n_doubs"])
    c_df.to_csv("/sc/arion/projects/CommonMind/pnm/genesis_U24/A6_Ruzicka/A6_Ruzicka_vir_demux_stats_wo_gt.tsv", sep='\t', index=False)


if __name__ == "__main__":
    main()
