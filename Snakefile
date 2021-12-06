import os
#from collections import OrderedDict
import glob2, re, math
import pandas as pd
from helper_py_scripts import update_logs
from snakemake.utils import validate
from itertools import repeat



###
#
# How to execute this script:
# ml anaconda3
# ml R/4.1.0
# conda activate snakemake
# ml star/2.7.5b
# snakemake --profile prachu_lsf --restart-times 2
#
###

###
# File Structure Example:
# fastq
# ── round2
# │   ├── Sample_NPSAD-20201013-A1-B-HTO
# │   │   └── fastq
# │   │       ├── checksum
# │   │       ├── NPSAD-20201013-A1-B-HTO_TCTCGCGC_HLM3WDSXY_L003_001.R1.fastq.gz
# │   │       └── NPSAD-20201013-A1-B-HTO_TCTCGCGC_HLM3WDSXY_L003_001.R2.fastq.gz
# │   ├── Sample_NPSAD-20201013-A1-cDNA
# │   │   └── fastq
# │   │       ├── checksum
# │   │       ├── NPSAD-20201013-A1-cDNA_CCGATGGTCT_HLM3WDSXY_L003_001.R1.fastq.gz
# │   │       ├── NPSAD-20201013-A1-cDNA_CCGATGGTCT_HLM3WDSXY_L003_001.R2.fastq.gz
#
###

# Create LSF profile as mentioned in https://github.com/Snakemake-Profiles/lsf
# The log files for the cluster execution are, by default, created in the structure (with respect to the dir of execution):
# error file: log/cluster/{rule name}/{wildcards}/{jobid}_{random_string}.err
# std log file: log/cluster/{rule name}/{wildcards}/{jobid}_{random_string}.out
#
# When encountering error for --limitSjdbInsertNsj in the STARsolo_sort rule read the associated Log.out, which is present as:
# /sc/arion/projects/psychAD/STARsolo_bams/{{num}}/Sample_{{id1}}-{{val}}/{{id1}}-{{id2}}_Log.out
# and the error line will be like:
# SOLUTION: re-run with at least --limitSjdbInsertNsj 1621424
# from which I extract 1621424
# grep -oP "^SOLUTION.*limitSjdbInsertNsj.* \K\w+" NPSAD-20201106-A2-AAGGGCCGCA-GAGGAATCAG_HVVKTDSXY_L001_001_Log.out
# module load everything before the start in the order:
#

configfile: "new_config.yaml"
#validate(config, "config.schema.json")  # Path to the scefic schema

#These are ordered
#round_num=[]
sample_name=[]
#sample_type=[]
#tag=[]


# Limitting Step for the run of Snakemake
with open(config['select_fastqs']) as fq:
    for line in fq:
        line_sp = line.split('/')
        #round_num.append(line_sp[0])
        #fn_sp = line_sp[0].split('-')
        #tag.append(line_sp[2].replace(fn_sp[0]+'_', '').strip())
        sample_name.append(line_sp[0].replace('-cDNA', ''))


# This df is to execute the demultiplexing part wherein we divide each set into its donor(SubID)
#df = pd.read_csv('/sc/arion/projects/psychAD/upload2synpase/final_NPSAD_snRNAseq_metrics_combo.csv')
#df.drop_duplicates('cDNA_ID', inplace=True)

#ruleorder: STARsolo_sort > run_calico_solo

#localrules: all, update_stats

# This variable controls the structure of the outputs from the rules produced for each set
# "" for no structure
# fold_struct="{id1}-cDNA/{id1}-cDNA" #"{num}/Sample_{id1}-cDNA/{id1}-cDNA"

# This variable controls the structure of the outputs from kallisto-bustools pipeline (rules create_FB, create_mismatch_fasta, build_kallisto_index, run_kallisto,
# run_bustools_correct, run_bustools_sort, run_bustools_count and also includes the rules cellSNP and vireoSNP) produced for each set
# # "" for no structure
# fold_struct_kb="{id1}-HTO/" #"{num}/Sample_{id1}-cDNA/"

# This variable controls the structure of the outputs from the rules create_h5ad_bustools, run_calico_solo and demux_samples_MULTIseq_solo_STARsolo produced for each set
# "" for no structure
# fold_struct_demux="{id1}-HTO/{id1}-HTO" #"{num}_Sample-{id1}-cDNA"



# Select final rule from the DAG
# Independent, standalone rules available: STARsolo* (STARsolo+filter+index; soon to be one), STARsolo*+PICARD tools (either or both of RNAseq metrics 
# and GC bias metrics), STARsolo*+kallisto_bustools+demux or all (STARsolo*+kb+calico_solo+PICARD/one/both). Dependent on the set of rules to be executed, update_stats will produce logs

# Value for the key "last_step" in config.json to run:
# 1) the whole pipeline - "all"
# 2) just STARsolo* - "STARsolo"
# 3) STARsolo* + PICARD's RNAseq metrics - "STARsolo_rnaseqmet"
# 4) STARsolo* + PICARD's GC bias metrics - "STARsolo_gcbiasmet"
# 5) STARsolo* + kallisto_bustools + demux by calico_solo/hashsolo - "STARsolo_kb_solo"
# 6) STARsolo* + PICARD's RNAseq metrics + PICARD's GC bias metrics - "STARsolo_PICARD"
# 7) cellSNP + vireoSNP - "STARsolo_gt_demux" (Not yet implemented)


# Check conditions to produce an aggregate 'log' of all measures
def check_log_version(conf_f) -> "list":
    map_file = conf_f['meta_data']
    mf_df = pd.read_csv(map_file, delimiter="\t", names=conf_f['meta_data_headers'])
    log_file = conf_f['log_all_stats']

    # If log file is not present then execute the rule to produce it
    if not os.path.isfile(log_file):
        # print("The file {} is not present!".format(log_file))
        return log_file

    else:
        lf_df = pd.read_csv(log_file, sep = "\t", header=[0, 1, 2])


    # The function "check_latest_columns" returns "True" if the columns are determined to be the same
    # in the latest version of the logs_file and the last version of the logs_file (this is the "map_file")
    # The extra "2" is for the extra columns added after reading all log files i.e. "doublet percent" and "negative percent"
    if update_logs.get_latest_extra_columns() + mf_df.shape[0] + 2 != lf_df.shape[1]:
        # print("The file {} has fewer columns than expected!\nRemoving the file and producing a newer version of the file".format(log_file))
        os.remove(log_file)
        return log_file

    # This function returns "True" if the last version of the logs_file (this is the "map_file") has all samples
    # present in the config['select_fastqs']
    if not all(sample for sample in sample_name if sample in lf_df["LAB"]["SAMPLE"]["SAMPLE"]):
        # print("The file {} doesn't contain all samples present in the \"fastq_files.txt\"!\nRemoving the file and producing a newer version of the file".format(log_file))
        os.remove(log_file)
        return log_file

    return []




# Run STARsolo_sort, index_bams and filter_GeneFull_STARsolo
def targets_STARsolo(conf_f) -> "str":

    inp_pref=conf_f['bams_dir']
    bai_suff=conf_f['bai']
    gf_mat=conf_f['genefull_lun_matrix']
    gf_features=conf_f['genefull_lun_features']
    gf_barcodes=conf_f['genefull_lun_barcodes']
    folder_st=conf_f['fold_struct']

    target_list = [f"{inp_pref}{folder_st}{bai_suff}", f"{inp_pref}{folder_st}{gf_mat}", f"{inp_pref}{folder_st}{gf_features}", f"{inp_pref}{folder_st}{gf_barcodes}"]
    return target_list



def targets_PICARD(conf_f, progs='all') -> "list":
    
    inp_pref = conf_f['bams_dir']
    rna_seq_suff = conf_f['rnaseq_metrics']
    gc_met_suff = conf_f['gc_bias_metrics']
    gc_summ_suff = conf_f['gc_summary_metrics']
    folder_st=conf_f['fold_struct']
    # bai_suff = conf_f['bai']
    # gf_mat = conf_f['genefull_lun_matrix']
    # gf_features = conf_f['genefull_lun_features']
    # gf_barcodes = conf_f['genefull_lun_barcodes']

    files_dict = {'RNAseq':rna_seq_suff, 'GC':[gc_met_suff, gc_summ_suff]}
    target_list = targets_STARsolo(conf_f=conf_f)
    if progs == 'all':
       # Run both PICARD programs
        for k, val in files_dict.items():
            if isinstance(val, list):
                for l_val in val:
                    target_list.append(f"{inp_pref}{folder_st}{l_val}")

            else:
                target_list.append(f"{inp_pref}{folder_st}{val}")

    elif progs == "RNAseq":
        # Only RNAseq metrics
        target_list.extend( [f"{inp_pref}{folder_st}{val}"  for val in files_dict[progs]] )

    elif progs == "GC":
        # Only GC bias metrics
        target_list.extend( [f"{inp_pref}{folder_st}{val}"  for val in files_dict[progs]] )
    else:
        # print("Wrong input; Check for editing errors!")
        return []

    return target_list



# To run STARsolo* + kb pipeline + (optional)PICARD progs

def targets_all(conf_f, PICARD=True, progs='all') -> "list":

    demuxed_mat_dir = conf_f['final_count_matrix_dir']
    demuxed_info_dir = conf_f['demultiplex_info_dir']    
    demux_mat_suff = conf_f['final_count_matrix_h5ad']
    demux_info_suff = conf_f['demultiplex_info']
    folder_st_bam=conf_f['fold_struct']
    folder_st=conf_f['fold_struct_demux']

    # STARsolo* + PICARD (any) progs
    if PICARD and (progs == 'all' or progs == 'RNAseq' or progs == 'GC'):
        target_list = targets_PICARD(conf_f=conf_f, progs=progs)

    else: # For any wrong values to progs or progs == None or just PICARD == False or everything else, run STARsolo*
        target_list = targets_STARsolo(conf_f=conf_f)


    # kb pipeline
    target_list.extend( [f"{demuxed_mat_dir}{folder_st}{demux_mat_suff}", f"{demuxed_info_dir}{folder_st}{demux_info_suff}"] )
    
    return target_list


# Add wildcards info here in the 'expand' function
def produce_targets(conf_f) -> "list":

    target_step = conf_f['last_step']
    folder_st_bam=conf_f['fold_struct']
    folder_st=conf_f['fold_struct_demux']


    if target_step == "all":
        target_files = targets_all(conf_f=conf_f, PICARD=True, progs=target_step)

        final_target_list = [expand(f"{target}", id1=sample_name) for target in target_files]

    elif target_step == "STARsolo_kb_solo":
        target_files = targets_all(conf_f=conf_f, PICARD=False, progs=None)

        final_target_list= [expand(f"{target}", id1=sample_name) for target in target_files]

    elif target_step == "STARsolo":
        target_files = targets_STARsolo(conf_f=conf_f)

        final_target_list= [expand(f"{target}", id1=sample_name) for target in target_files]
        
    elif target_step == "STARsolo_PICARD":
        target_list = targets_PICARD(conf_f=conf_f, progs='all')

        final_target_list= [expand(f"{target}", id1=sample_name) for target in target_files]

    elif target_step == "STARsolo_rnaseqmet":
        target_list = targets_PICARD(conf_f=conf_f, progs='RNAseq')

        final_target_list= [expand(f"{target}", id1=sample_name) for target in target_files]

    elif target_step == "STARsolo_gcbiasmet":
        target_list = targets_PICARD(conf_f=conf_f, progs='GC')

        final_target_list= [expand(f"{target}", id1=sample_name) for target in target_files]

    # Not yet implemented
    elif target_step == "STARsolo_gt_demux":
        pass        

    else:
        # print("Wrong inputs to produce_targets function!!")
        return []
        
    final_target_list.append(check_log_version(conf_f=conf_f))

    return final_target_list
        

# Add wildcards info here in the 'expand' function
def stats_produce_inp(wildcards):
    STAR_log=expand(f"{config['bams_dir']}{config['fold_struct']}{config['STAR_log_final']}", id1=wildcards.id1)
    SS_G_Feat=expand(f"{config['bams_dir']}{config['fold_struct']}{config['gene_features']}", id1=wildcards.id1)
    SS_GF_Feat=expand(f"{config['bams_dir']}{config['fold_struct']}{config['genefull_features']}", id1=wildcards.id1)
    SS_G_Summ=expand(f"{config['bams_dir']}{config['fold_struct']}{config['gene_summary']}", id1=wildcards.id1)
    SS_GF_Summ=expand(f"{config['bams_dir']}{config['fold_struct']}{config['genefull_summary']}", id1=wildcards.id1)
    SS_Barcodes=expand(f"{config['bams_dir']}{config['fold_struct']}{config['barcodes_stats']}", id1=wildcards.id1)
    PICARD_GC=expand(f"{config['bams_dir']}{config['fold_struct']}{config['gc_summary_metrics']}", id1=wildcards.id1)
    PICARD_RNAseq=expand(f"{config['bams_dir']}{config['fold_struct']}{config['rnaseq_metrics']}", id1=wildcards.id1)
    Demultiplex_info=expand(f"{config['demultiplex_info_dir']}{config['fold_struct_demux']}{config['demultiplex_info']}", id1=wildcards.id1)

    if config['last_step'] == "all":
        
        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_GC, PICARD_RNAseq, Demultiplex_info]

    elif config['last_step'] == "STARsolo_kb_solo":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, Demultiplex_info]

    elif config['last_step'] == "STARsolo":
        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes]

    elif config['last_step'] == "STARsolo_PICARD":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_GC, PICARD_RNAseq]

    elif config['last_step'] == "STARsolo_rnaseqmet":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_RNAseq]

    elif config['last_step'] == "STARsolo_gcbiasmet":

        return [STAR_log, SS_G_Feat, SS_GF_Feat, SS_G_Summ, SS_GF_Summ, SS_Barcodes, PICARD_GC]

    # Not yet implemented
    elif config['last_step'] == "STARsolo_gt_demux":
        pass

    else:
        raise ValueError(f"Invalid Value to {config['last_step']} in the yaml file")



def stats_produce_params(wildcards, input):

    # Files for update log files
    log_dict = {'STAR_log':[config['STAR_log_final'], '--ss_l'], 'PICARD_GC':[config['gc_summary_metrics'], '--pc_gc'], 'PICARD_RNAseq':[config['rnaseq_metrics'], '--pc_rs'], 'SS_G_Feat':[config['gene_features'], '--ss_g_f'],
              'SS_GF_Feat':[config['genefull_features'], '--ss_gf_f'], 'SS_G_Summ':[config['gene_summary'], '--ss_g_s'], 'SS_GF_Summ':[config['genefull_summary'], '--ss_gf_s'], 'SS_Barcodes':[config['barcodes_stats'], '--ss_bc'],
               'Demultiplex_info':[config['demultiplex_info'], '--dem_info']}

    cons_param = ""
    for i in range(len(input)):
        s_temp = pd.Series(input[i])
        for k, v in log_dict:
            # All files should contain the pattern string
            if all(s_temp.str.contains(v[0])):
                curr_param = v[1] + f" {{input[{i}]}}"
                cons_param += curr_param


    return cons_param




# print(produce_targets(conf_f=config))

rule all:
    input:
        produce_targets(conf_f=config)
        # check_log_version()




# To check if the last digit of the line in the error log of STARsolo is a number
def check_isnumber(x):
    try:
        int(x)
        return True
   
    except ValueError:
        return False



# # This function reads the log file created per attempt to change the parameter "limitsjdbInsertNsj" in STARsolo   
# # We can use this to similarly change other params in the log
# def get_limitsjdbval(wildcards, resources):
#     # This is to check the log file produced after each attempt for the error value
#     file_p_temp = fold_struct.format(id1=wildcards.id1)
#     log_list = glob2.glob("{}{}_STARsolo_log.txt*".format(config['bams_dir'], file_p_temp))
#     for log_file in log_list: 
#         with open(log_file) as fin:
#             for line in fin:
#                 if line.startswith("SOLUTION") and "limitSjdbInsertNsj" in line and check_isnumber(line.split()[-1]):
#                     print("Found an Error with the parameter limitSjdbInsertNsj. Changing from the default value of 1000000 to {}".format(line.split()[-1]))
#                     return line.split()[-1]

#                 else:
#                     continue

#             # Empty but existing file
#             else:
#                 continue
                    

#    # This is to check the parameters file (if there was a previous successful run)
#     else:
#         if os.path.isfile("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)):
#             with open("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)) as fin:
#                 for line in fin:
#                     print("Found limitSjdbInsertNsj value from the previous successfull run in {}. Using the same value".format(config['star_params_dir']))
#                     return re.search("--limitSjdbInsertNsj ([0-9]+) ", line).group(1)

#         return 1000000


#     return 1000000


# def get_limitsjcollapsed(wildcards, resources):
#     # This is to check the log file produced after each attempt for the error value
#     file_p_temp = fold_struct.format(id1=wildcards.id1)
#     log_list = glob2.glob("{}{}_STARsolo_log.txt*".format(config['bams_dir'], file_p_temp))
#     print(log_list)
#     for log_file in log_list:
#         with open(log_file) as fin:
#             for line in fin:
#                 print(line)
#                 if line.startswith("SOLUTION") and "limitSjdbInsertNsj" in line and check_isnumber(line.split()[-1]):
#                     print("Found an Error with the parameter limitSjdbInsertNsj. Changing the value of limitOutSJcollapsed from the default value of 1000000 to {}".format(line.split()[-1]))
#                     return line.split()[-1]
           
#                 elif line.startswith("Solution") and "limitOutSJcollapsed" in line:
#                     print("Found an Error with limitOutSJcollapsed. Changing from the default value of 1000000 to {}".format(1000000*(1+resources.attempt)))
#                     return 1000000*(1+resources.attempt)

#                 else:
#                     print("Else block inside For block")
#                     continue

#             # Empty but existing file
#             else:
#                 continue


#    # This is to check the parameters file (if there was a previous successful run)
#     else:
#         if os.path.isfile("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)):
#             with open("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)) as fin:
#                 for line in fin:
#                     print("Found limitOutSJcollapsed value from the previous successfull run in {}. CHangin it to {}".format(config['star_params_dir'], 
#                         int(re.search("--limitOutSJcollapsed ([0-9]+) ", line).group(1))*(1+resources.attempt)))
#                     return int(re.search("--limitOutSJcollapsed ([0-9]+) ", line).group(1))*(1+resources.attempt)

#         print("Else block compliment to For block")
#         return 1000000


    # return 1000000




# We can use this to similarly change other params in the log
def get_limitsjdbval_coll(wildcards, resources):
    '''
    This function reads the log file created per attempt to change the parameter "limitsjdbInsertNsj" and "limitOutSJcollapsed" in STARsolo
    Serially produce output as a list in the sequence "limitsjdbInsertNsj", "limitOutSJcollapsed", etc.
    '''
    # This is to check the log file produced after each attempt for the error value
    file_p_temp = config['fold_struct'].format(id1=wildcards.id1)
    log_list = glob2.glob("{}{}_STARsolo_log.txt*".format(config['bams_dir'], file_p_temp))
    ins_nsj = 1000000
    sj_collap = 1000000
    for log_file in log_list: 
        with open(log_file) as fin:
            for line in fin:
                if line.startswith("SOLUTION") and "limitSjdbInsertNsj" in line and check_isnumber(line.split()[-1]):
                    # print("Found an Error with the parameter limitSjdbInsertNsj. Changing defaulkt values of the parameters \
                        # \"limitSjdbInsertNsj\" and \"limitOutSJcollapsed\" from the default value of 1000000 to {}".format(line.split()[-1]))
                    ins_nsj = line.split()[-1] if line.split()[-1] > ins_nsj else ins_nsj
                    sj_collap = ins_nsj

                elif line.startswith("Solution") and "limitOutSJcollapsed" in line:
                    # print("Found an Error with limitOutSJcollapsed. Changing from the default value of 1000000 to {}".format(1000000*(1+resources.attempt)))
                    sj_collap = 1000000*(1+resources.attempt)
                    ins_nsj = sj_collap

                else:
                    continue

            # Empty but existing file
            else:
                continue
                    

   # This is to check the parameters file (if there was a previous successful run)
    else:
        if os.path.isfile("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)):
            with open("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)) as fin:
                for line in fin:
                    # print("Found values of \"limitSjdbInsertNsj\" and \"limitOutSJcollapsed\" from the previous successfull run in {}. Using the same value".format(config['star_params_dir']))
                    ins_nsj = re.search("--limitSjdbInsertNsj ([0-9]+) ", line).group(1)
                    sj_collap = re.search("--limitOutSJcollapsed ([0-9]+) ", line).group(1)

        # return 1000000


    # return 1000000
    return [ins_nsj, sj_collap]




def get_log_file(wildcards, resources):

    file_p_temp = config['fold_struct'].format(id1=wildcards.id1)
    return f"{config['bams_dir']}{file_p_temp}_STARsolo_log.txt_{resources.attempt}"


# def get_hto_demux(wildcards):
#     hto_d_f = [f for f in glob2.glob("/sc/arion/projects/psychAD/Single_cell_data/preprocess/{a}/HTO_demultiplexed/Sample_{b}_HTOdemux.csv".format(a=wildcards.num, b=wildcards.id1.replace('-cDNA', '-*HTO*')))]
#     return hto_d_f



# def get_multiseq(wildcards):
#     multi_f = [f for f in glob2.glob("/sc/arion/projects/psychAD/Single_cell_data/preprocess/{a}/HTO_demultiplexed/Sample_{b}_MULTIseq.csv".format(a=wildcards.num, b=wildcards.id1.replace('-cDNA', '-*HTO*')))]
#     return multi_f



def get_filt_barcodes(wildcards):
    barc_f = [f for f in glob2.glob("{}{}{a}_{b}.txt".format(config['filt_barcodes_dir'], config['filt_barcodes'], a=wildcards.num.replace('round_num', ''), b=wildcards.id1))]
    return barc_f



# def get_filt_barcodes_cr(wildcards):
#     barc_f = [f for f in glob2.glob("/sc/arion/projects/psychAD/pnm/inp_for_vireo/filt_bc_vireo_{a}_{b}_cr.txt".format(a=wildcards.num.replace('round_num', ''), b=wildcards.id1))]
#     return barc_f



def calc_donors(wildcards):
    temp_df = pd.read_csv(config['wet_lab_info'])
    n_donors = temp_df.loc[temp_df[config['columns_to_pick'][0]] == wildcards.id1, config['columns_to_pick'][2]].shape[0]
    return n_donors



def get_donors(wildcards):
    temp_df = pd.read_csv(config['meta_data_geno_samp'], skiprows=1, names=["SubID", "Orig_VCF_ID", "Samples"])
    donors = ','.join(temp_df.loc[temp_df["Samples"] == wildcards.id1[:-6], "Orig_VCF_ID"].to_list())
    return donors



# Mem allocation per thread (here thread is core/tasks), everything is in MB
def allocate_mem_PM(wildcards, attempt):
    return 2500*attempt+2500


def allocate_mem_KBP(wildcards, attempt):
    return 2500*attempt+2500


def allocate_mem_DXP(wildcards, attempt):
   return 3500*attempt+3500


def allocate_mem_cS(wildcards, attempt):
   return 1000*attempt+1500


def allocate_mem_vS(wildcards, attempt):
   return 10000*attempt+10000



rule STARsolo_sort:
    input:
        R1=f"{config['cDNA_fastqs_dir']}{config['fold_struct']}{config['R1_suffix']}",
        R2=f"{config['cDNA_fastqs_dir']}{config['fold_struct']}{config['R2_suffix']}"

    priority: 10

    params:
        gtf=config['gtf_file'],
        genome_dir=config['genome_dir'],
        #output_prefix=lambda wildcards, input: input[0][:-29],
        overhang=config['sjdboverhang'],
        opt_params=get_limitsjdbval_coll,
        # limitsjdbval=get_limitsjdbval,
        chemistry=config['soloType'], # For STARsolo
        whitelist=config['whitelist'], # V3 whitelist
        UMI_length=config['umi_len'], # V3 
        SAM_attr=config['SAM_attr'],
        features=config['features'],
        save_params=f"{config['star_params_dir']}Sample_{{id1}}-cDNA.txt",
        star_def_log_out=f"{config['bams_dir']}{config['fold_struct']}_Log.out",
        # limitsjcollap=get_limitsjcollapsed,
        solo_cell_filter=config['solo_cell_filter'],
        out_pref=lambda wildcards, output: output[7][:-13]
        #cell_filtering not present in 2.7.5b
        #cell_filtering="EmptyDrops_CR" # Cell Filtering matching CellRanger 3.0.0 from Lun et. al. 2019
        # Few more parameters to match CellRanger >=4.0.0

    output:
        f"{config['bams_dir']}{config['fold_struct']}{config['bam']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['STAR_log_final']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['genefull_features']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['genefull_summary']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['gene_features']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['gene_summary']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['barcodes_stats']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['genefull_lun_matrix']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['genefull_lun_features']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['genefull_lun_barcodes']}"


    log:
        f"{config['bams_dir']}{config['fold_struct']}_STARsolo_log.txt"
   
    resources:
        mem_mb=70000,
        time_min=1440,
        attempt=lambda wildcards, attempt: attempt

    threads: 7

    group: "BAM_processing"

    run:
        shell(
        """
        ml {config[STAR_version]}
        echo "{params.opt_params[0]}, {params.opt_params[1]}, {resources.attempt}"
        if [ ! -d {config[star_params_dir]} ]; then mkdir -p {config[star_params_dir]}; fi
        STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.opt_params[0]} --twopassMode Basic --readFilesCommand zcat --readFilesIn {input.R2} {input.R1} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} --limitOutSJcollapsed {params.opt_params[1]} --outSAMtype BAM SortedByCoordinate --runThreadN 12 --outFileNamePrefix {config[bams_dir]}{params.struct_fold}_ &> {log}_{resources.attempt}
        gzip {params.out_pref}*
        a=$(grep -n "^##### Final effective command line" {params.star_def_log_out} | cut -d ":" -f1)
        a=$((a+1))
        if [ ! -f "{params.save_params}" ]; then 
               tail -n +${{a}} {params.star_def_log_out} | head -1 > {params.save_params}
        else 
               tail -n +${{a}} {params.star_def_log_out} | head -1 > {params.save_params}_{resources.attempt}
               cmp --silent {params.save_params} {params.save_params}_{resources.attempt} && rm {params.save_params}_{resources.attempt} || rm {params.save_params} && mv {params.save_params}_{resources.attempt} {params.save_params}
        fi
        """
        )



rule index_bams:
    input:
        f"{config['bams_dir']}{config['fold_struct']}{config['bam']}"

    priority: 10

    output:
        f"{config['bams_dir']}{config['fold_struct']}{config['bai']}"

    threads: 1

    resources:
        mem_mb=10000,
        time_min=1440

    group: "BAM_processing"

    shell:
        """
        ml samtools
        samtools index {input}
        """



rule Picard_GC_bias_metrics:
    input:
        bams=f"{config['bams_dir']}{config['fold_struct']}{config['bam']}"

    priority: 9

    output:
        f"{config['bams_dir']}{config['fold_struct']}{config['gc_bias_metrics']}",
        f"{config['bams_dir']}{config['fold_struct']}{config['gc_summary_metrics']}"

    params:
        output_pref=lambda wildcards, output: output[0].replace(f"{config['gc_bias_metrics']}", '_'),
        window_size=config["window_size"],
        genome_fasta=config["genome_fasta"]

    resources:
        mem_mb=allocate_mem_PM,
        time_min=120
       
    threads: 2

    group: "PICARD_metrics"
     
    shell:
        """
        ml picard
        java -jar $PICARD CollectGcBiasMetrics I={input.bams} O={output[0]} CHART={params.output_pref}gc_bias_metrics.pdf SCAN_WINDOW_SIZE={params.window_size} S={params.output_pref}summary_metrics.txt R={params.genome_fasta}
        """
      


rule Picard_RNAseq_metrics:
    input:
        bams=f"{config['bams_dir']}{config['fold_struct']}{config['bam']}"

    priority: 9

    output:
        f"{config['bams_dir']}{config['fold_struct']}{config['rnaseq_metrics']}"

    params:
        flat_ref=config['flat_ref'],
        strand=config['strand']

    resources:
        mem_mb=allocate_mem_PM,
        time_min=420

    threads: 4

    group: "PICARD_metrics"

    shell:
        """
        ml picard
        java -jar $PICARD CollectRnaSeqMetrics I={input.bams} O={output} REF_FLAT={params.flat_ref} STRAND={params.strand}
        """



rule create_FB:
    input:
        R1="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=config['fold_struct'].replace('-cDNA', '-HTO'), suff=config['R1_suffix']),
        R2="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=config['fold_struct'].replace('-cDNA', '-HTO'), suff=config['R2_suffix'])

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['feature_barcodes']}"

    resources:
        mem_mb=100, #allocate_mem_KBP,
        time_min=10

    shell: 
        """
        python3 helper_py_scripts/create_Feat_Barc.py
        """

         


rule create_mismatch_fasta:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['feature_barcodes']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['features_mismatch_fa']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['features_mismatch_t2g']}"

    params:
        headers=config['headers'] # Does the feature barcodes file hasve headers

    resources:
        mem_mb=100, #allocate_mem_KBP,
        time_min=10


    shell: 
        """
        if [[ "{params.headers}" == "yes" ]]; then
            python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]} --header

        else
            python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]}

        fi
        """


rule build_kallisto_index:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['features_mismatch_fa']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['features_mismatch_idx']}"

    resources:
        mem_mb=2000, #allocate_mem_KBP,
        time_min=10

    shell:
        """
        ml kallisto
        kallisto index -i {output} -k 15 {input}
        """


rule run_kallisto:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['features_mismatch_idx']}",
        R1="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=config['fold_struct'].replace('-cDNA', '-HTO'), suff=config['R1_suffix']),
        R2="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=config['fold_struct'].replace('-cDNA', '-HTO'), suff=config['R2_suffix'])

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_file']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['ec_matrix']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['tx']}"

    params:
        output_pref= lambda wildcards, output: output[0].replace(f"{config['bus_file']}", ''),
        chemistry=config['chemistry']

    threads: 1

    resources:
        mem_mb=lambda wildcards, attempt: 1000*attempt + 500, #allocate_mem_KBP,
        time_min=40

    shell:
        """
        ml kallisto
        kallisto bus -i {input[0]} -o {params.output_pref} -x {params.chemistry} -t 8 {input.R1} {input.R2}
        """


rule run_bustools_correct:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_file']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_file_corrected']}"

    params:
        whitelist=config['whitelist']

    resources:
        mem_mb=3000, #allocate_mem_KBP,
        time_min=10

    shell:
        """
        ml bustools
        bustools correct -w {params.whitelist} {input} -o {output}
        """

        
rule run_bustools_sort:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_file_corrected']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_file_sorted']}"

    threads: 1

    resources:
        mem_mb=6000, #allocate_mem_KBP,
        time_min=10

    shell:
        """
        ml bustools
        bustools sort -t 4 -o {output} {input}
        """


rule run_bustools_count:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_file_sorted']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['ec_matrix']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['tx']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['features_mismatch_t2g']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_count_dir']}{config['bus_count_mtx']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_count_dir']}{config['bus_count_genes']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_count_dir']}{config['bus_count_barcodes']}"

    params:
        output_pref=lambda wildcards, output: output[0].replace(f"{config['bus_count_mtx']}", '')

    resources:
        mem_mb=3000, #allocate_mem_KBP,
        time_min=10

    shell:
        """
        ml bustools
        bustools count -o {params.output_pref} --genecounts -g {input[3]} -e {input[1]} -t {input[2]} {input[0]}
        """


rule create_h5ad_bustools:
    input:
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_count_dir']}{config['bus_count_mtx']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_count_dir']}{config['bus_count_genes']}",
        f"{config['kallisto_bustools_dir']}{config['fold_struct_kb']}{config['bus_count_dir']}{config['bus_count_barcodes']}"

    priority: 10

    output:
        f"{config['h5ad_bustools_dir']}{config['fold_struct_demux']}.h5ad"

    resources:
        mem_mb=3000, #allocate_mem_KBP,
        time_min=10

    shell:
        """
        python3 helper_py_scripts/create_h5ad_from_bustools.py {input[0]} {input[1]} {input[2]} -o {output}
        """



rule run_calico_solo:
    input:
        f"{config['h5ad_bustools_dir']}{config['fold_struct_demux']}.h5ad",
        starsolo_out=f"{config['bams_dir']}{config['fold_struct']}{config['genefull_lun_matrix']}" #get_STARsolo_mat

    priority: 8
   
    output:
        f"{config['calico_solo_dir']}{config['fold_struct_demux']}{config['calico_solo_h5ad']}"
  
    params:
        mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        min_genes=config['min_genes_per_cell'], # Min #genes per cell
        min_cells=config['min_cells_per_gene']  # Min #cells expressing a gene for it to pass the filter

    threads: 2

    resources:
        mem_mb=5000,
        time_min=10

    shell: 
        """
        python3 helper_py_scripts/create_h5ad_from_calico_solo.py {input[0]} {input[1]} {output} -m {params.mito} -g {params.min_genes} -c {params.min_cells}
        sleep 100
        """



rule demux_samples_MULTIseq_solo_STARsolo:
    input:
        f"{config['calico_solo_dir']}{config['fold_struct_demux']}{config['calico_solo_h5ad']}",
        starsolo_out=f"{config['bams_dir']}{config['fold_struct']}{config['genefull_lun_matrix']}"

    priority: 8

    output:
        f"{config['final_count_matrix_dir']}{config['fold_struct_demux']}{config['final_count_matrix_h5ad']}",
        f"{config['demultiplex_info_dir']}{config['fold_struct_demux']}{config['demultiplex_info']}"

    params:
        mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        min_genes=config['min_genes_per_cell'], # Min #genes per cell
        min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
        cols=config['columns_to_pick'],  # Which columns of the wet lab info file correspond RESPECTIVELY to cDNA_ID(should correspond to the name of the processed file), HTO numbers and Donors/SubIDs (Header names and not numbers)
        genes_info=config['gene_info_file'] # File containing gene names and gene ids for annotations

    
    resources:
        mem_mb=allocate_mem_DXP,
        time_min=30

    shell:
        """
        python3 helper_py_scripts/demul_samples.py {input[0]} {input[1]} {output[0]} {output[1]} {params.genes_info} {params.samples_info} -m {params.mito} -g {params.min_genes} -c {params.min_cells} --headers {params.cols}
        sleep 100
        """
        
       


# UMI tag is turned on. Therefore, PCR duplicates are included
rule cellSNP:
    input:
        bc=get_filt_barcodes,
        bams=f"{config['bams_dir']}{config['fold_struct']}{config['bam']}"

    output:
        f"{config['cellsnp_dir']}{config['fold_struct_kb']}{config['cellsnp_cells']}",
        f"{config['cellsnp_dir']}{config['fold_struct_kb']}{config['cellsnp_base']}"

    params:
        ref_snps=config['ref_snps'],
        umi_tag=config['umi_tag'],
        cell_tag=config['cell_tag'],
        output_prefix=lambda wildcards, output: output[0].replace(f"/{config['cellsnp_cells']}", '')

    threads: 8

    resources:
        mem_mb=allocate_mem_cS

    shell: "cellsnp-lite -s {input.bams} -b {input.bc} -O {params.output_prefix} -R {params.ref_snps} -p {config[n_proc]} --minMAF {config[min_maf]} --minCOUNT {config[min_aggr_count]} --cellTAG {params.cell_tag} --UMItag {params.umi_tag} --genotype --gzip"


rule vireoSNP:
    input:
        f"{config['cellsnp_dir']}{config['fold_struct_kb']}{config['cellsnp_cells']}",
        config['genotyped_vcf'],
        config['meta_data_geno_samp']
  
    output:
        f"{config['vireosnp_dir']}{config['fold_struct_kb']}{config['donors_vcf']}",
        f"{config['filt_vcf_dir']}{config['fold_struct_kb']}{config['filt_vcf']}"

    params:
        n_donors=calc_donors,
        donors=get_donors,
        geno_tag=config['donor_genotype'],
        output_prefix=lambda wildcards, output: output[0].replace(f"/{config['donors_vcf']}", '')

    threads: 7

    resources:
        mem_mb=allocate_mem_vS
        
    shell:
        """
        ml bcftools
        bcftools view {input[1]} -R {input[0]} -s {params.donors} -Oz -o {output[1]}
        vireo -c {input[0]} -d {output[1]} -N {params.n_donors} -o {params.output_prefix} -t {params.geno_tag} --noPlot

        """



#rule cellSNP_cr:
   #input:
    #    "/sc/arion/projects/psychAD/Single_cell_data/alignment/{num}/cDNA/results/Sample_{id1}/outs/possorted_genome_bam.bam",
#        bc=get_filt_barcodes_cr

 #  output:
  #      "/sc/arion/projects/psychAD/demux_gt/cellSNP/{num}/Sample_{id1}_cr/{config['cellsnp_cells']}",
  #      "/sc/arion/projects/psychAD/demux_gt/cellSNP/{num}/Sample_{id1}_cr/{config['cellsnp_base']}"

  # params:
   #     ref_snps="/sc/arion/projects/psychAD/Single_cell_data/cellSNP_ref/CMC_SNParray_cellSNP_ref.vcf.gz",
    #    umi_tag="UB",
     #   cell_tag="CB",
      #  output_prefix=lambda wildcards, output: output[0].replace(f"/{config['cellsnp_cells']}", '')

  # threads: 8

   #resources:
    #    mem_mb=allocate_mem_cS

   #shell: "cellsnp-lite -s {input[0]} -b {input.bc} -O {params.output_prefix} -R {params.ref_snps} -p 20 --minMAF 0.1 --minCOUNT 20 --cellTAG {params.cell_tag} --UMItag {params.umi_tag} --genotype --gzip"



#rule vireoSNP_cr:
 #  input:
 #       "/sc/arion/projects/psychAD/demux_gt/cellSNP/{num}/Sample_{id1}_cr/{config['cellsnp_cells']}",
 #       config['genotyped_vcf'],
 #       config['meta_data_geno_samp']
#
#   output:
#        "/sc/arion/projects/psychAD/demux_gt/vireoSNP/{num}/Sample_{id1}_cr/{config['donors_vcf']}",
#        "/sc/arion/projects/psychAD/demux_gt/bcf_filt_vcf/{num}/Sample_{id1}_cr/{config['filt_vcf']}"
#
 #  params:
  #      n_donors=calc_donors,
  #      donors=get_donors,
  #      geno_tag="GT",
  #      output_prefix=lambda wildcards, output: output[0].replace(f"/{config['donors_vcf']}", '')
#
 #  threads: 7

  # resources:
   #     mem_mb=allocate_mem_vS

   #shell:
   #     """
   #     ml bcftools
   #     bcftools view {input[1]} -R {input[0]} -s {params.donors} -Oz -o {output[1]}
   #     vireo -c {input[0]} -d {output[1]} -N {params.n_donors} -o {params.output_prefix} -t {params.geno_tag} --noPlot
#
 #       """

#rule genetic_check:
   #input:



rule split_bams:
    input:
        bam=f"{config['bams_dir']}{config['fold_struct']}{config['bam']}",
        barcode_list="" # Barcodes vs sample name txt file, produced after producing final_count_matrices

    params:
        output_pref=f"{config['split_bams_dir']}{config['fold_struct']}_"
    output:
        f"{config['split_bams_proxy_dir']}{{num}}_Sample_{{id1}}.txt" # Proxy to the output

    shell:
        """
        helper_py_scripts/split_bam_indiv_barc_samtools.sh {input.bam} {input.barcode_list} {params.output_pref} > {output}
        sleep 60
        """



# Rule to run update.py
rule update_stats:
    input:
        inp_files = stats_produce_inp
            
    params:
        map_file=config['meta_data'],
        cl_inp_params=stats_produce_params

    output:
        config['log_all_stats']

    threads: 1

    resources:
        mem_mb=500,
        time_min=30

    shell:
        """
        python3 helper_py_scripts/update_logs.py {params.cl_inp_params} -m {params.map_file} -o {output}
        sleep 60
        """


#rule 
