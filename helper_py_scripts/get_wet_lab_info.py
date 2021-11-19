#!/usr/bin/env python

import pandas as pd
import re, os
from openpyxl import load_workbook
import subprocess
import sys



fin=sys.argv[1]

wb2 = pd.read_csv('/sc/arion/projects/psychAD/upload2synpase/final_NPSAD_snRNAseq_metrics_combo.csv')

if fin == 'All_HTOs_info.csv':

   try:
       pd.read_csv('/sc/arion/projects/psychAD/pnm/All_HTOs_info.csv', sep='\t', header=False)
   
   except: 
       col_vals = pd.unique(wb2[['hashtag', 'ab..barcode']].values.ravel('K'))
       a = np.split(cols_vals, 2)
       ht_num = [ v.replace('#','') for v in a[0].tolist() ]
       ht_seq = a[1].tolist()
       hto_df = pd.DataFrame({'HTO_name':ht_num, 'HTO_seq':ht_seq})
       hto_df["HTO_name"] = hto_df["HTO_name"].map(lambda x: 'HTO'+str(x))
       hto_df.to_csv('/sc/arion/projects/psychAD/pnm/All_HTOs_info.csv', sep='\t', header=False, index=False)


elif fin == 'All_Samples_HTO_info.csv':

   try:
       pd.read_csv('/sc/arion/projects/psychAD/pnm/All_Samples_HTO_info.csv', sep='\t', header=False)

   except:
       col_vals = pd.unique(wb2[['set_ID', 'hashtag']].values.ravel('K'))
       a = np.split(cols_vals, 2)
       ht_num = [ v.replace('#','') for v in a[1].tolist() ]
       samp_names = [ v.replace("_", "-") for v in a[0].tolist() ]
       sample_hto_df = pd.DataFrame({'Sample':samp_names, 'HTOs':ht_num})
       sample_hto_df["HTOs"] = sample_hto_df["HTOs"].map(lambda x: 'HTO'+str(x))
       sample_hto_df.to_csv('/sc/arion/projects/psychAD/pnm/All_Samples_HTO_info.csv', sep='\t', header=False, index=False)

else:
   raise ValueError("Listed file cannot be created by this script!")

