#!/usr/bin/env python3
#import sys

"""Create Feature Barcodes file for KITEseq.

This script is used to create featurebarcodes csv file given a file 
containing sample names, its corresponding hto names and hto barcodes, respectively.
"""

# Examples
# --------
#     $ python3 create_Feat_Barc.py samples_info.txt -o sample_1_fb.csv -s sample_1 -c set_name hto_name hto_bc

# samples_info.txt

# set_name hto_name hto_bc
# sample_1 hto1 ATCTATGGTTG
# sample_1 hto3 ATGAATGGTTG
# sample_2 hto3 ATCCGGTGTTG
# sample_1 hto4 ATCTATGGTTG
# sample_2 hto2 AGCTATGGTTG
# sample_2 hto6 ATTTATGGTTG
# sample_1 hto5 ATCTCTTGTTG
# sample_2 hto8 ATCTACTATTG
# ...
# ...
# ...
# ...

# sample_1_fb.csv
# sample_1 hto1 ATCTATGGTTG
# sample_1 hto3 ATGAATGGTTG
# sample_1 hto5 ATCTCTTGTTG
# sample_1 hto4 ATCTATGGTTG

# Help
# -----
#     $ python3 create_Feat_Barc.py -h

# Usage
# ------
#     $ python3 create_Feat_Barc.py <input_file> -o <output_file> -s <sample_name> -c <space-separated list of 3 columns>


import pandas as pd
import argparse, os
from time import sleep


def get_argument_parser():
    """Generate and return argument parser."""
    #Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Create FeatureBarcodes "
        "file for a given sample from a wet lab manifest"
        )

    # Positional Paramaters
    parser.add_argument('donor_file', help="Path to file that contains "
        "HTO info for each set (wet_lab_file created by create_wet_lab_info.py file)"
        )

    # Optional parameters
    parser.add_argument('-o', '--output', help="Name of the output file (csv file)")
    parser.add_argument('-s', '--sample_name', help="Sample name "
        "(should be present in the wet lab file)")
    parser.add_argument('-c', '--columns', nargs=3, 
        help="List of column names RESPECTIVELY to pool_name(should "
        "contain the sample name provided to this script), HTO names "
        "and corresponding HTO sequence.", metavar=('pool_name', 'HTO_name', 'HTO_barcode'), 
        default=['unique_sample_ID', 'hashtag', 'ab_barcode'],
        )

    return parser


def main(): 
    """Main function"""
    
    parser = get_argument_parser()
    args = parser.parse_args()

    # Store output_file and Create necessary folders
    op = args.output

    # Create parent dir(s) to the output
    if not os.path.isdir(op.replace('/' + os.path.basename(op), '')):
        os.makedirs(op.replace(os.path.basename(op), ''))
                
    fout = op if op.endswith('.csv') else op + '.csv'
    # print(fout)

    #val_list = []
    samp = args.sample_name
    cols = args.columns


    pd_df = pd.read_csv(args.donor_file, sep='\t')
    a = pd_df.loc[pd_df[cols[0]] == samp, cols[1:]]
    # a[cols[1]] = a[cols[1]].apply(lambda x: x.replace('HTO#', ''))

    # Custom-requirement
    # If htos are present in one-line
    # example:
    # 1,2,3 ATGCTAGCTAG,ATCGATGCTG,GATCGATCGT
    if a.shape[0] > 1:
        htos=a.iloc[0,0].split(',')
        hto_barc_seq=a.iloc[0,1].split(',')

        df_to_save=[]
        for i in range(len(htos)):
            df_to_save.append( ('HTO' + htos[i], hto_barc_seq[i]) )

        df = pd.DataFrame(df_to_save)
    
        df.to_csv(fout, sep=',', header=False, index=False, mode='w+')
    else:
        a.to_csv(fout, sep=',', header=False, index=False, mode='w+')


if __name__ == '__main__':
    main()
    sleep(10)