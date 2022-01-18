#!/usr/bin/env python
#import sys
import pandas as pd
import re, argparse, os

#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Create FeatureBarcodes file for a given sample from a wet lab manifest")

# Positional Paramaters
parser.add_argument('wet_lab_file', help="Path to file that contains HTO info for each set")

# Optional parameters
parser.add_argument('-o', '--output', help="Name of the output file")
parser.add_argument('-s', '--sample_name', help="Name of the sample")
parser.add_argument('-c', '--columns', nargs=3, help="List of column names RESPECTIVELY to set_ID(should contain the sample name provided to this script), \
    HTO numbers and Donors/SubIDs.", metavar=('set_ID', 'HTO_name', 'HTO_barcode'), default=['set_ID', 'hto', 'hto_barcode'])


args = parser.parse_args()

fout = args.output
print(fout)
#samples_hto_df = pd.read_csv(snakemake.input[1], sep='\t', names=['Samples', 'HTOs'], header=None)
#all_htos = pd.read_csv(snakemake.input[2], sep='\t', names=['HTO_name', 'HTO_seq'], header=None)
#samples_hto_df['HTOs'] = samples_hto_df['HTOs'].apply(eval)

#val_list = []
samp_val = args.sample_name.replace('-', '_') +'_HTO'
cols = args.columns


#for htos in samples_hto_df[samp_val == samples_hto_df['Samples']]['HTOs']:
 #  for hto in sorted([int(x) for x in htos]):
  #     seq = all_htos[all_htos['HTO_name'] == 'HTO' + str(hto)]['HTO_seq'].values[0]
   #    seq = seq + 'A' if len(seq)%2 == 0 else seq
    #   val_list.append(('HTO' + str(hto), seq))


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
