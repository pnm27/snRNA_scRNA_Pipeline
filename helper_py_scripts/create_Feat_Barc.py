#!/usr/bin/env python
#import sys
import pandas as pd
import re

fout = snakemake.output[0]
print(fout)
#samples_hto_df = pd.read_csv(snakemake.input[1], sep='\t', names=['Samples', 'HTOs'], header=None)
#all_htos = pd.read_csv(snakemake.input[2], sep='\t', names=['HTO_name', 'HTO_seq'], header=None)
#samples_hto_df['HTOs'] = samples_hto_df['HTOs'].apply(eval)

#val_list = []
samp_val = re.search('\/Sample_(NPSAD.*)-cDNA\/', fout).group(1)

#for htos in samples_hto_df[samp_val == samples_hto_df['Samples']]['HTOs']:
 #  for hto in sorted([int(x) for x in htos]):
  #     seq = all_htos[all_htos['HTO_name'] == 'HTO' + str(hto)]['HTO_seq'].values[0]
   #    seq = seq + 'A' if len(seq)%2 == 0 else seq
    #   val_list.append(('HTO' + str(hto), seq))


pd_df = pd.read_csv("/sc/arion/projects/psychAD/upload2synpase/final_NPSAD_snRNAseq_metrics_combo.csv")
a = pd_df.loc[pd_df['set_ID'] == samp_val].loc[:, [ 'hashtag', 'ab..barcode']].reset_index(drop=True)
a.columns = ['HTO', 'HTO_seq']
a['HTO'] = a.HTO.map(lambda x: 'HTO'+x.replace('#', ''))
a.to_csv(fout, sep=',', header=False, index=False, mode='w+')
