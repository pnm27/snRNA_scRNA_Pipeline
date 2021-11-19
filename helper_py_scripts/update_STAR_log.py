#!/usr/bin/env python

import pandas as pd
import os
import itertools
import numpy as np
import glob2
import re


def write_logs(log_file, big_df, mapper):
   temp_df = pd.read_csv(log_file, names=["cols", "vals"], delimiter=r"|", skiprows=[7, 22, 27, 34])
   temp_df["vals"] = temp_df.vals.str.strip()
   temp_df["cols"] = temp_df.cols.str.strip()
   col_keys = [['STAR'], mapper.curr_val.tolist()]
   new_row = []
   new_row.append(re.search('/round([0-9]+)/', log_file).group(1))
   new_row.append(os.path.basename(log_file).replace("_Log.final.out", ""))
   # Add values to a list in the same sequence as the final output file/dataframe
   for prog, val in big_df.columns.tolist():
       if prog != "Sample" and prog != "Batch":
           add_value = temp_df.loc[temp_df["cols"] == mapper.loc[mapper["curr_val"] == val, "val_in_log"].values[0], "vals"].values[0]
           new_row.append(add_value.replace(" ","/"))
   
   return new_row
   

# Iterate on all log file
all_logs = glob2.glob("/sc/arion/projects/psychAD/STARsolo_bams/round*/*/*_Log.final.out")
#print(len(all_logs))
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_gc_bias_metrics.txt", sep='\t', skiprows=range(6))
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_rnaseq_metrics.txt", sep='\t', nrows=1, skiprows=6)
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_summary_metrics.txt", sep='\t', skiprows=6)
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_qual_score_dist.txt", sep='\t', skiprows=7)


map_names = pd.read_csv("Final_out_MAP.csv", delimiter="\t", names=["val_in_log", "curr_val"])
for i in range(len(all_logs)):
   cols = [["STAR"], map_names.curr_val.tolist()]
   cl = list(itertools.product(*cols))
   cl = [('Batch', 'Round'), ('Batch', 'Preparer'), ('Sample', 'Sample')]+cl
   try:
       #if os.path.isfile(snakemake.output[0]) :
       combo_log = pd.read_csv(snakemake.output[0], sep = "\t", header=[0, 1])
       #if not(combo_log["Sample"]["Sample"].str.contains(os.path.basename(all_logs[i]).replace("_Log.final.out", "")).any()) :
           #combo_log.loc[len(combo_log.iloc[:, 0])] = write_logs(all_logs[i], combo_log, map_names)
   
       #combo_log.to_csv(snakemake.output[0], sep = "\t", index=False)
   except:
       #with open(snakemake.output[0], 'w+') as fout:
       combo_log = pd.DataFrame(columns=pd.MultiIndex.from_tuples(cl, names=["prog", "value"]))
       #combo_log.loc[0] = write_logs(all_logs[i], combo, map_names)
       #df.to_csv(snakemake.output[0], sep = "\t", index=False)

   finally:
       #print(type(combo_log))
       #print(combo_log.columns)
       if not(combo_log["Sample"]["Sample"].str.contains(os.path.basename(all_logs[i]).replace("_Log.final.out", "")).any()) :
           combo_log.loc[len(combo_log.iloc[:, 0])] = write_logs(all_logs[i], combo_log, map_names)

       combo_log.to_csv(snakemake.output[0], sep = "\t", index=False)

