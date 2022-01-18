#!/usr/bin/env python

import pandas as pd, numpy as np
import os, re, glob2, argparse

def ret_samp(x):
    x = x[1:]
    
    return x



def ret_samp_names(y, df_info):
    #print(len(y), len(bc))
    if df_info.loc[df_info['SubID'] == y, 'cDNA_ID'].isnull().all():
        return ""
    else:
        return '-'.join(df_info.loc[df_info['SubID'] == y, 'cDNA_ID'].to_list()[0].split('-')[:3])[:-1]



#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Produce inputs (text files) per sample for cellSNP, containing filtered barcodes for each")

parser.add_argument('geno_info', help="Genotype info file (Expected to be csv formatted)")
parser.add_argument('sample_meta_data', help="Sample meta data (Expected to be csv formatted)")
parser.add_argument('geno_vcf', help="")

# Optional parameters
parser.add_argument('-o', '--output', help="Name of the output file.")


args = parser.parse_args()

temp_df = pd.read_csv('VA_BR_Clinical_combo_matched_CMC_AMPAD_10082020.csv')
temp_df = temp_df.loc[temp_df["SNP_report:Genotyping_Sample_ID"].dropna().index.tolist(), :]
temp_df.reset_index(drop=True)

test_run = temp_df.loc[temp_df['SubNum'].isin(values) , ["SubNum", "SNP_report:Genotyping_Sample_ID"]]
test_run['SubNum'] = test_run['SubNum'].map(lambda x: 'M'+str(x))
test_run.reset_index(drop=True, inplace=True)

vcf = pd.read_csv('CMC_MSSM_NPSAD_genotyped_all.vcf.gz', skiprows=46, compression='gzip', sep='\t')


df_shan = pd.read_csv('/sc/arion/projects/psychAD/upload2synpase/final_NPSAD_snRNAseq_metrics_combo.csv')

values = list(map(ret_samp, df_shan['SubID'].tolist()))
values = [ int(v) for v in values ]


headers = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + test_run["SNP_report:Genotyping_Sample_ID"].tolist()

headers = [h for h in headers if h in vcf.columns]
headers = [h for i, h in enumerate(headers) if h not in headers[:i]]

vcf_sub = vcf.loc[:, headers]

test_run['Samp_ID'] = test_run['SubNum'].apply(ret_samp_names, args=(df_shan, ))
test_run.drop_duplicates(inplace=True, ignore_index=True)

test_run.to_csv('Hash_map_gtID_SubID.csv', index=False)
