#!/usr/bin/env python
#import sys

"""Create Feature Barcodes file for KITEseq.

This script is used to create featurebarcodes csv file given a file 
containing sample names, its corresponding hto names and hto barcodes, respectively.

Help
-----
    $ python3 create_Feat_Barc.py -h

Usage
------
    $ python3 create_Feat_Barc.py <input_file> -o <output_file> -s <sample_name> -c <space-separated list of 3 columns>

Examples
--------
    $ python3 create_Feat_Barc.py samples_info.txt -o sample_1_fb.csv -s sample_1 -c set_name hto_name hto_bc

samples_info.txt

set_name hto_name hto_bc
sample_1 hto1 ATCTATGGTTG
sample_1 hto3 ATGAATGGTTG
sample_2 hto3 ATCCGGTGTTG
sample_1 hto4 ATCTATGGTTG
sample_2 hto2 AGCTATGGTTG
sample_2 hto6 ATTTATGGTTG
sample_1 hto5 ATCTCTTGTTG
sample_2 hto8 ATCTACTATTG
...
...
...
...

sample_1_fb.csv
sample_1 hto1 ATCTATGGTTG
sample_1 hto3 ATGAATGGTTG
sample_1 hto5 ATCTCTTGTTG
sample_1 hto4 ATCTATGGTTG

"""
import pandas as pd
import re, argparse, os

#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Create FeatureBarcodes file for a given sample from a wet lab manifest")

# Positional Paramaters
parser.add_argument('wet_lab_file', help="Path to file that contains HTO info for each set")

# Optional parameters
parser.add_argument('-o', '--output', help="Name of the output file (csv file)")
parser.add_argument('-s', '--sample_name', help="Name of the sample")
parser.add_argument('-c', '--columns', nargs=3, help="List of column names RESPECTIVELY to set_ID(should contain the sample name provided to this script), \
    HTO numbers and Donors/SubIDs.", metavar=('set_ID', 'HTO_name', 'HTO_barcode'), default=['set_ID', 'hto', 'hto_barcode'])


args = parser.parse_args()

fout = args.output if args.output.endswith('.csv') else args.output + '.csv'
print(fout)

#val_list = []
samp_val = args.sample_name.replace('-', '_') +'_HTO'
cols = args.columns


pd_df = pd.read_csv(args.wet_lab_file, sep='\t')
a = pd_df.loc[pd_df[cols[0]] == samp_val, cols[1:]] # Series as an output, HTOs at a[0] and barcode sequences at a[1]
a[cols[1]] = a[cols[1]].apply(lambda x: x.replace('HTO#', ''))

# Custom-requirement
htos=a.iloc[0,0].split(',')
hto_barc_seq=a.iloc[0,1].split(',')

df_to_save=[]
for i in range(len(htos)):
    df_to_save.append( ('HTO' + htos[i], hto_barc_seq[i]) )

df = pd.DataFrame(df_to_save)
df.to_csv(fout, sep=',', header=False, index=False, mode='w+')
