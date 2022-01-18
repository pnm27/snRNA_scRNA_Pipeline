#!/usr/bin/env python

import pandas as pd
import os
import itertools
import numpy as np
import glob2
import re
from time import sleep
import errno


# This function is to read files with extension '.stats', which are formatted weirdly
# and return a pandas Dataframe for easy use
def get_df(inp_path):
   col1 = []
   col2 = []
   with open(inp_path) as f1:
       for line in f1:
           col1.append(line.strip().split()[0])
           col2.append(line.strip().split()[1])
   
   n_df = pd.DataFrame({'cols':col1,'vals':col2})
   return n_df


# This function is used to caluclate ratios like doublet pct and negative pct
def calc_ratio(numer, denom):
    ratios = []
    for i in range(len(numer)):
        ratios.append(int(numer[i])/int(denom[i]))
        
    return ratios


# Main function to write rows of samples into the
def write_logs(big_df, mapper, all_files_dict, **kwargs):
   
   new_row = []
   new_row.append(kwargs["round_num"])
   new_row.append(kwargs["sample"])
   new_row.append(kwargs["set_num"])
   new_row.append(kwargs["prep"])
   new_row.append(kwargs["rep"])
   # Store Barcode.stats and Feature.stats files as DF in a list for easy access
   # Seq is Feature.stats for Gene, Feature.stats for GeneFull, Barcode.stats
   stats_file = [x for k, x in all_files_dict.items() if '.stats' in x ]
   list_df = [get_df(x) for x in stats_file ]
   
   # Add values to a list in the same sequence as the final output file/dataframe
   for prog, sub_prog, val in big_df.columns.tolist():
       if prog != "LAB":
           if prog == "STAR":
               temp_df = pd.read_csv(all_files_dict["STAR_final"], names=["cols", "vals"], delimiter=r"|", skiprows=[7, 22, 27, 34])
               temp_df["vals"] = temp_df.vals.str.strip()
               temp_df["cols"] = temp_df.cols.str.strip()
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[mapper["curr_val"] == val, "val_in_log"].values[0], "vals"].values[0]
               new_row.append(add_value.replace(" ","/"))

           elif sub_prog == "GC":
               temp_df = pd.read_csv(all_files_dict["PICARD_GC"], sep='\t', skiprows=6)
               add_value = temp_df.loc[0, mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GC"), "val_in_log"].values[0]]
               new_row.append(add_value)
               
           elif sub_prog == "RNASEQMETRIC":
               temp_df = pd.read_csv(all_files_dict["PICARD_RNASeq"], sep='\t', nrows=1, skiprows=6)
               add_value = temp_df.loc[0, mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "RNASEQMETRIC"), "val_in_log"].values[0]]
               new_row.append(add_value)

           elif sub_prog == "GENE_FEATURE":
               temp_df = get_df(all_files_dict["Gene_Features"]) 
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENE_FEATURE"), "val_in_log"].values[0], "vals"].values[0]
               new_row.append(add_value)

           elif sub_prog == "GENE_SUMM":
               temp_df = pd.read_csv(all_files_dict["Gene_Summary"], names=['cols', 'vals'])
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENE_SUMM"), "val_in_log"].values[0], "vals"].values[0]
               new_row.append(add_value)

           elif sub_prog == "GENEFULL_FEATURE":
               temp_df = get_df(all_files_dict["GeneFull_Features"])
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENEFULL_FEATURE"), "val_in_log"].values[0], "vals"].values[0]
               new_row.append(add_value)

           elif sub_prog == "GENEFULL_SUMM":
               temp_df = pd.read_csv(all_files_dict["GeneFull_Summary"], names=['cols', 'vals'])
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENEFULL_SUMM"), "val_in_log"].values[0], "vals"].values[0]
               new_row.append(add_value)

           elif sub_prog == "BARCODE_STATS":
               temp_df = get_df(all_files_dict["Barcodes_stats"])
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "BARCODE_STATS"), "val_in_log"].values[0], "vals"].values[0]
               new_row.append(add_value)
           
           elif sub_prog == "DEMUX":
               temp_df = pd.read_csv(all_files_dict["Demultiplex_stats"], names=['cols', 'vals'], skiprows=1, sep='\t')
               add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "DEMUX"), "val_in_log"].values[0], "vals"].values[0]
               if val == "N_CELLS_AFTER_DEMUX_CS":
                   add_value = add_value[:add_value.find('Name:')]
                   add_value= re.sub('[^\S\r\n]+', ':', add_value)
                   add_value = add_value[:-1]
                   new_row.append(add_value)

               else:
                   new_row.append(add_value)           

           else:
               raise ValueError(f'This extra column exists in the output file-All_logs.csv: {prog}, {sub_prog}, {val}')


   return new_row


# Mention the file structure for each log file
#'/sc/arion/projects/psychAD/STARsolo_bams/round2/Sample_NPSAD-20201117-C2-cDNA/NPSAD-20201117-C2-cDNA_L004_001_Solo.out/Gene/Features.stats'
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_rnaseq_metrics.txt", sep='\t', nrows=1, skiprows=6)
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_summary_metrics.txt", sep='\t', skiprows=6)
#pd.read_csv("../STARsolo_bams/round2/Sample_NPSAD-20201125-A1-cDNA/NPSAD-20201125-A1-cDNA_qual_score_dist.txt", sep='\t', skiprows=7)
#demultiplex/solo/round1_Sample-NPSAD-20201105-A1-cDNA_calicosolo_info.tsv
   


# Extract the round info and Sample info from fastqs: val_set is a list of sets, which comprises of (round, sample_name)

# Constraints to limit the fastq output
# no_op = [('1', 'NPSAD-20201021-A1-cDNA', 'A', '1'), ('3', 'NPSAD-20201218-A1-cDNA', 'A', '1'), ('1', 'NPSAD-20201013-A1-cDNA', 'A', '1'), ('3', 'NPSAD-20201218-A2-cDNA', 'A', '2'),  \
#          ('1', 'NPSAD-20201016-A1-cDNA', 'A', '1'), ('1', 'NPSAD-20201022-A1-cDNA', 'A', '1'), ('1', 'NPSAD-20201021-A2-cDNA', 'A', '2'), ('1', 'NPSAD-20201022-A2-cDNA', 'A', '2')]
# val_set = [(re.search('/round([0-9]+)/', f).group(1), re.search('/Sample_(NPSAD-.*-cDNA)/', f).group(1), os.path.basename(f).split('-')[2][0], os.path.basename(f).split('-')[2][1]) for f in glob2.glob("/sc/arion/projects/psychAD/fastqs/round*/*-cDNA/*R1.fastq.gz") if 'round4' not in f]
# val_set = [vals for vals in val_set if vals not in no_op]


#all_logs = glob2.glob("/sc/arion/projects/psychAD/STARsolo_bams/round*/*/*_Log.final.out")
#print(len(all_logs))

new_cols_to_add = [['ROUND', 'LAB', 'BATCH'], ['SAMPLE', 'LAB', 'SAMPLE'], ['SET', 'LAB', 'BATCH'], ['PREPARER', 'LAB', 'BATCH'], ['REP', 'LAB', 'BATCH']]

# Function to conditionally run this script through Snakemake if the current file has fewer columns that the last version of this script
def get_latest_extra_columns():
    global new_cols_to_add
    # if curr_length == len(new_cols_to_add)
    #     return True
    # else
    #     return False
    return len(new_cols_to_add)


# Run this only when executed through Snakemake
if __name__ == "__main__":
    print("This script will update info/logs for all Samples simulatneously!!!")
    print("Run this script after all logs described below has been produced.")
    print("Logs produced by STARsolo-Log.final.out, Summary.csv and Features.stats for Gene and GeneFull, info from PICARD's CollectGCBiasMetrics, RNASeqMetrics and demultiplex output stats for each sample will be used.")
    print("If any of the above mentioned files (9) does not exist that sample won't be processed")
    #print("Currently all files exist for only for samples in round1-3 excluding:")
    #print("round1_Sample-NPSAD-20201021-A2\nround3_Sample-NPSAD-20201218-A1\nround1_Sample-NPSAD-20201013-A1\nround1_Sample-NPSAD-20201022-A1\n\
    #round3_Sample-NPSAD-20201218-A2\nround1_Sample-NPSAD-20201016-A1\nround1_Sample-NPSAD-20201022-A2\nround1_Sample-NPSAD-20201021-A1")
    print("\n\n")
    print("Naming Convention Assumed to run this script")
    print("Fastq files: NPSAD-<yyyymmdd>-<one_letter_for_preparer><one_digit_for_rep_num>-cDNA_<Lane_info>.R<1/2>.fastq.gz")
    print("\t\t\tExample: NPSAD-20201110-A1-cDNA_L001_001.R1.fastq.gz")

    out=snakemake.output[0]
    map_names = pd.read_csv(snakemake.params["map_file"], delimiter="\t", names=["val_in_log", "curr_val", "prog", "sub_prog", "desc"])
    #cl = pd.DataFrame([['Round', 'lab', 'Batch'], ['Preparer', 'lab', 'Batch'], ['Sample', 'lab', 'Sample']], columns=list(map_names.columns.values)[1:-1])
    #cols = map_names.iloc[:, 1:-1]
    #cols = cols.append(cl, ignore_index=True)
    #cols = cols.iloc[:, [0, 2, 1]]
    cl = pd.DataFrame(new_cols_to_add, columns=list(map_names.columns.values)[1:-1])
    cols = map_names.iloc[:, 1:-1]
    cl = cl.append(cols, ignore_index=True)

    # Change column header order to: curr_val, sub_prog, prog so that the output file looks like:
    #Prog
    #Sub_prog
    #curr_val

    #eg:
    #STARsolo
    #GENE_SUMM
    #N_READS
    cl = cl.iloc[:, [1, 2, 0]]

    # Integrity check
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["PICARD_GC"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["PICARD_RNAseq"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["SS_G_Feat"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["SS_GF_Feat"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["SS_G_Summ"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["SS_GF_Summ"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["SS_Barcodes"])
    assert len(snakemake.input["STAR_log"]) == len(snakemake.input["Demultiplex_info"])

    # If file doesn't exist create one else open as pandas dataframe
    try:
        #if os.path.isfile(snakemake.output[0]) :
        combo_log = pd.read_csv(out, sep = "\t", header=[0, 1, 2])

        # Catch older files that may have lesser columns than expected
        if combo_log.shape[1] != cl.shape[0]:
            combo_log = pd.DataFrame(columns=pd.MultiIndex.from_frame(cl, names=["prog", "sub_prog", "curr_val"]))


    except:
        #with open(snakemake.output[0], 'w+') as fout:
        combo_log = pd.DataFrame(columns=pd.MultiIndex.from_frame(cl, names=["prog", "sub_prog", "curr_val"]))


    for i in range(len(snakemake.input["STAR_log"])):
        

           # flag
        flag = 1
        #cols = [["STAR"], map_names.curr_val.tolist()]
        #cl = list(itertools.product(*cols))
        #cl = [('lab', 'Batch', 'Round'), ('lab', 'Batch', 'Preparer'), ('Sample', 'Sample', 'Sample')]+cl
           #if os.path.isfile(snakemake.output[0]) :
           #if not(combo_log["Sample"]["Sample"].str.contains(os.path.basename(all_logs[i]).replace("_Log.final.out", "")).any()) :
               #combo_log.loc[len(combo_log.iloc[:, 0])] = write_logs(all_logs[i], combo_log, map_names)

           #combo_log.to_csv(snakemake.output[0], sep = "\t", index=False)
           #with open(snakemake.output[0], 'w+') as fout:
           #combo_log = pd.DataFrame(columns=pd.MultiIndex.from_tuples(cl, names=["prog", "sub_prog", "value"]))
           #combo_log.loc[0] = write_logs(all_logs[i], combo, map_names)
           #df.to_csv(snakemake.output[0], sep = "\t", index=False)

        #print(type(combo_log))
        #print(combo_log.columns)


        # files_list = ["/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_Log.final.out".format(num=val_set[i][0], samp_name=val_set[i][1]), 
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_summary_metrics.txt".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_rnaseq_metrics.txt".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_Solo.out/Gene/Features.stats".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_Solo.out/GeneFull/Features.stats".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_Solo.out/Gene/Summary.csv".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_Solo.out/GeneFull/Summary.csv".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/STARsolo_bams/round{num}/Sample_{samp_name}/{samp_name}_Solo.out/Barcodes.stats".format(num=val_set[i][0], samp_name=val_set[i][1]),
        #               "/sc/arion/projects/psychAD/demultiplex/info/round{num}_Sample-{samp_name}_STARsolo_info.tsv".format(num=val_set[i][0], samp_name=val_set[i][1])]

        files_dict = {"STAR_final":snakemake.input["STAR_log"][i], "PICARD_GC":snakemake.input["PICARD_GC"][i], "PICARD_RNASeq":snakemake.input["PICARD_RNAseq"][i], "Gene_Features":snakemake.input["SS_G_Feat"][i], "GeneFull_Features":snakemake.input["SS_GF_Feat"][i], 
                      "Gene_Summary":snakemake.input["SS_G_Summ"][i], "GeneFull_Summary":snakemake.input["SS_GF_Summ"][i], "Barcodes_stats":snakemake.input["SS_Barcodes"][i], "Demultiplex_stats":snakemake.input["Demultiplex_info"][i]}

        # file_check = list(map(os.path.isfile, files_list))

        sample_name = re.search('/(NPSAD-.*)_Log.final.out', files_dict["STAR_final"]).group(1)
        r_num = int(re.search('/round([0-9]+)/', files_dict["STAR_final"]).group(1))
        preparer = sample_name.split('-')[2][0]
        replicate = sample_name.split('-')[2][1]
        set_val = sample_name[:-6]

        # Redundant file check
        # if not all(map(os.path.isfile, files_list)):
        #     for key, val in files_dict.items():
        #         raise FileNotFoundError(errno.ENOENT, "From round {} the sample {} does not have the {} file at: ".format(round_num, sample_name, key), val)


           
        if not(combo_log['LAB']['SAMPLE']['SAMPLE'].str.contains(sample_name).any()) :
           #combo_log.loc[len(combo_log.iloc[:, 0])] = write_logs(combo_log, map_names, file1, file2, file3, file4, file5, file6, file7, file8, round_num=val_set[i][0], sample=val_set[i][1], prep=val_set[i][2], rep=val_set[i][3])
            combo_log.loc[len(combo_log.iloc[:, 0])] = write_logs(combo_log, map_names, files_dict, round_num=r_num, sample=sample_name, prep=preparer, rep=replicate, set_num=set_val)

        print(f"Finished adding {sample_name} to the file")

    try:
        combo_log["STARsolo", "DEMUX", "DOUBLET_PCT"] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_DOUBLET_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
    except:
        print("Doublet ratio for: {}".format(sample_name))

    try:
        combo_log["STARsolo", "DEMUX", "NEGATIVE_PCT"] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_NEGATIVE_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
    except:
        print("Negative ratio for: {}".format(sample_name))


    combo_log.to_csv(out, sep = "\t", index=False)


    sleep(30)
