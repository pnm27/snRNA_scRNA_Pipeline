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

configfile: "config.json"
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
fold_struct="{id1}-cDNA/{id1}-cDNA" #"{num}/Sample_{id1}-cDNA/{id1}-cDNA"

# This variable controls the structure of the outputs from kallisto-bustools pipeline (rules create_FB, create_mismatch_fasta, build_kallisto_index, run_kallisto,
# run_bustools_correct, run_bustools_sort, run_bustools_count and also includes the rules cellSNP and vireoSNP) produced for each set
# "" for no structure
fold_struct_kb="{id1}-cDNA/" #"{num}/Sample_{id1}-cDNA/"

# This variable controls the structure of the outputs from the rules create_h5ad_bustools, run_calico_solo and demux_samples_MULTIseq_solo_STARsolo produced for each set
# "" for no structure
fold_struct_demux="{id1}-cDNA/{id1}-cDNA"#"{num}_Sample-{id1}-cDNA"



# Select final rule from the DAG
# Independent, standalone rules available: STARsolo* (STARsolo+filter+index; soon to be one), STARsolo*+PICARD tools (either or both of RNAseq metrics 
# and GC bias metrics), STARsolo*+kallisto_bustools+demux or all (STARsolo*+kb+calico_solo+PICARD/one/both). Dependent on the set of rules to be executed, update_stats will produce logs

# Value for the key "last_step" in config.json to run:
# 1) the whole pipeline - "all"
# 2) just STARsolo* - "STARsolo"
# 3) STARsolo* + PICARD's RNAseq metrics - "STARsolo_rnaseqmet"
# 4) STARsolo* + PICARD's GC bias metrics - "STARsolo_gcbiasmet"
# 5) STARsolo* + kallisto_bustools + demux by calico_solo/hashsolo - "STARsolo_kb_solo"
# 6) cellSNP + vireoSNP (Not yet implemented)


# Run STARsolo_sort, index_bams and filter_GeneFull_STARsolo
def targets_STARsolo(inp_pref=config['bams_dir'], folder_st=fold_struct, bai_suff=config['bai'], gf_mat=config['genefull_lun_matrix'], gf_features=config['genefull_lun_features'], 
    gf_barcodes=config['genefull_lun_barcodes']) -> "str":

        # expand( zip,  id1=sample_name),
        # expand(f"{config['bams_dir']}{fold_struct}{}", zip,  id1=sample_name),
        # expand(f"{config['bams_dir']}{fold_struct}{}", zip,  id1=sample_name),
    target_list = [f"{inp_pref}{fold_struct}{bai_suff}", f"{inp_pref}{fold_struct}{gf_mat}", f"{inp_pref}{fold_struct}{gf_features}", f"{inp_pref}{fold_struct}{gf_barcodes}"]
    return target_list



def targets_PICARD(inp_pref=config['bams_dir'], progs='all',folder_st=fold_struct, rna_seq_suff=config['rnaseq_metrics'], gc_met_suff=config['gc_bias_metrics'], gc_summ_suff=config['gc_summary_metrics'], bai_suff=config['bai']
    , gf_mat=config['genefull_lun_matrix'], gf_features=config['genefull_lun_features'], gf_barcodes=config['genefull_lun_barcodes']) -> "list":

    files_dict = {'RNAseq':rna_seq_suff, 'GC':[gc_met_suff, gc_summ_suff]}
    target_list = targets_STARsolo(inp_pref=inp_pref, folder_st=folder_st, bai_suff=bai_suff, gf_mat=config['genefull_lun_matrix'], gf_features=config['genefull_lun_features'], gf_barcodes=config['genefull_lun_barcodes'])
    if progs == 'all':
       # Run both PICARD programs
        for k, val in files_dict.items():
            if isinstance(val, list):
                for l_val in val:
                    target_list.extend(f"{inp_pref}{l_val}")

            else:
                target_list.extend(f"{inp_pref}{val}")

    elif progs == "RNAseq":
        # Only RNAseq metrics
        target_list.extend(f"{inp_pref}{val}"  for val in files_dict[prog] )

    elif progs == "GC":
        # Only GC bias metrics
        target_list.extend(f"{inp_pref}{val}"  for val in files_dict[prog] )
    else:
        print("Wrong input; Check for editing errors!")
        return []

    return target_list



# To run STARsolo* + kb pipeline + (optional)PICARD progs

def targets_all(bam_pref=config['bams_dir'], bai_suff=config['bai'],folder_st_bam=fold_struct, demuxed_mat_dir=config['final_count_matrix_dir'], demuxed_info_dir=config['demultiplex_info_dir'], folder_st=fold_struct_demux, demux_mat_suff=config['final_count_matrix_h5ad']
    , demux_info_suff=config['demultiplex_info'], PICARD=True, progs='all', rna_seq_suff=config['rnaseq_metrics'], gc_met_suff=config['gc_bias_metrics'], gc_summ_suff=config['gc_summary_metrics'], gf_mat=config['genefull_lun_matrix'] 
    , gf_features=config['genefull_lun_features'], gf_barcodes=config['genefull_lun_barcodes']) -> "list":

    # STARsolo* + PICARD (any) progs
    if PICARD and ( progs == 'all' or progs == 'RNAseq' or progs == 'GC'):
        target_list = targets_PICARD(inp_pref=bam_pref, progs=progs, bai_suff=bai_suff, folder_st=folder_st_bam, rna_seq_suff=rna_seq_suff, gc_met_suff=gc_met_suff, gc_summ_suff=gc_summ_suff, gf_mat=config['genefull_lun_matrix']
            , gf_features=config['genefull_lun_features'], gf_barcodes=config['genefull_lun_barcodes'])

    else: # For any wrong values to progs or progs == None or just PICARD == False or everything else, run STARsolo*
        target_list = targets_STARsolo(inp_pref=bam_pref, folder_st=folder_st_bam, bai_suff=bai_suff, gf_mat=config['genefull_lun_matrix'], gf_features=config['genefull_lun_features'], gf_barcodes=config['genefull_lun_barcodes'])


    # kb pipeline
    target_list.extend([f"{demuxed_mat_dir}{folder_st}{demux_mat_suff}", f"{demuxed_info_dir}{folder_st}{demux_info_suff}"])
    print(f"{demuxed_mat_dir}{folder_st}{demux_mat_suff}")
    return target_list



#def targets_all(): -> "list"
    



def produce_targets(config_file=config, folder_st_bam=fold_struct, folder_st=fold_struct_demux, id=sample_name, PICARD_progs='all') -> "list":

    target_step = config_file['last_step']

    if target_step == "all" or target_step == "STARsolo_rnaseqmet" or target_step == "STARsolo_gcbiasmet":
        target_files = targets_all(bam_pref=config_file['bams_dir'], bai_suff=config_file['bai'],folder_st_bam=folder_st_bam, demuxed_mat_dir=config_file['final_count_matrix_dir'], demuxed_info_dir=config_file['demultiplex_info_dir'],
            folder_st=folder_st, demux_mat_suff=config_file['final_count_matrix_h5ad'], demux_info_suff=config_file['demultiplex_info'], PICARD=True, progs=PICARD_progs, rna_seq_suff=config_file['rnaseq_metrics'],
            gc_met_suff=config_file['gc_bias_metrics'], gc_summ_suff=config_file['gc_summary_metrics'], gf_mat=config_file['genefull_lun_matrix'], gf_features=config_file['genefull_lun_features'], gf_barcodes=config_file['genefull_lun_barcodes'])

        return [expand(f"{target}", id1=id) for target in target_files]

    if target_step == "STARsolo_kb_solo":
        target_files = targets_all(bam_pref=config_file['bams_dir'], bai_suff=config_file['bai'],folder_st_bam=folder_st_bam, demuxed_mat_dir=config_file['final_count_matrix_dir'], demuxed_info_dir=config_file['demultiplex_info_dir'],
            folder_st=folder_st, demux_mat_suff=config_file['final_count_matrix_h5ad'], demux_info_suff=config_file['demultiplex_info'], PICARD=False, progs=None, rna_seq_suff=config_file['rnaseq_metrics'],
            gc_met_suff=config_file['gc_bias_metrics'], gc_summ_suff=config_file['gc_summary_metrics'], gf_mat=config_file['genefull_lun_matrix'], gf_features=config_file['genefull_lun_features'], gf_barcodes=config_file['genefull_lun_barcodes'])

        return [expand(f"{target}", id1=id) for target in target_files]

    elif target_step == "STARsolo":
        target_files = targets_STARsolo(inp_pref=config_file['bams_dir'], folder_st=folder_st_bam, bai_suff=config_file['bai'], gf_mat=config_file['genefull_lun_matrix'], gf_features=config_file['genefull_lun_features'], 
            gf_barcodes=config_file['genefull_lun_barcodes'])

        return [expand(f"{target}", id1=id) for target in target_files]
        

    # Not yet implemented
    elif target_step == "cellSNP + vireoSNP":
        pass        

    else:
        print("Wrong inputs to produce_targets function!!", )
        return []
        


def check_log_version() -> "list":
    map_file = config['meta_data']
    mf_df = pd.read_csv(map_file, delimiter="\t", names=config['meta_data_headers'])
    log_file = config['log_all_stats']

    # If log file is not present then execute the rule to produce it
    if not os.path.isfile(log_file):
        print("The file {} is not present!".format(log_file))
        return log_file

    else:
        lf_df = pd.read_csv(log_file, sep = "\t", header=[0, 1, 2])


    # The function "check_latest_columns" returns "True" if the columns are determined to be the same
    # in the latest version of the logs_file and the last version of the logs_file (this is the "map_file")
    # The extra "2" is for the extra columns added after reading all log files i.e. "doublet percent" and "negative percent"
    if update_logs.get_latest_extra_columns() + mf_df.shape[0] + 2 != lf_df.shape[1]:
        print("The file {} has fewer columns than expected!\nRemoving the file and producing a newer version of the file".format(log_file))
        os.remove(log_file)
        return log_file

    # This function returns "True" if the last version of the logs_file (this is the "map_file") has all samples
    # present in the config['select_fastqs']
    if not all(sample for sample in sample_name if sample in lf_df["LAB"]["SAMPLE"]["SAMPLE"]):
        print("The file {} doesn't contain all samples present in the \"fastq_files.txt\"!\nRemoving the file and producing a newer version of the file".format(log_file))
        os.remove(log_file)
        return log_file

    return []


print(produce_targets(config_file=config, id=sample_name))

rule all:
    input:
        produce_targets(config_file=config, id=sample_name)
        # expand(f"{config['bams_dir']}{fold_struct}{config['bai']}", zip,  id1=sample_name),
        #expand(f"{config['bams_dir']}{{num}}/Sample_{{id1}}/{{id1}}{config['gene_lun_matrix']}", zip,  id1=sample_name),
        #expand(f"{config['bams_dir']}{{num}}/Sample_{{id1}}/{{id1}}{config['gene_lun_features']}", zip,  id1=sample_name),
        #expand(f"{config['bams_dir']}{{num}}/Sample_{{id1}}/{{id1}}{config['gene_lun_barcodes']}", zip,  id1=sample_name),
        # expand(f"{config['bams_dir']}{fold_struct}{config['genefull_lun_matrix']}", zip,  id1=sample_name),
        # expand(f"{config['bams_dir']}{fold_struct}{config['genefull_lun_features']}", zip,  id1=sample_name),
        # expand(f"{config['bams_dir']}{fold_struct}{config['genefull_lun_barcodes']}", zip,  id1=sample_name),
        #expand(f"{config['bams_dir']}{fold_struct}{config['gc_bias_metrics']}", zip,  id1=sample_name),
        #expand(f"{config['bams_dir']}{fold_struct}{config['summary_metrics.txt']}", zip,  id1=sample_name),
        #expand(f"{config['bams_dir']}{fold_struct}{config['rnaseq_metrics']}", zip,  id1=sample_name),
        #expand(f"{config['final_count_matrix_dir']}{fold_struct_demux}{config['bus_file_sorted']}", zip,  id1=sample_name),
        #expand(f"{config['demultiplex_info_dir']}{fold_struct_demux}{config['demultiplex_info']}", zip,  id1=sample_name)
        #expand(f"{config['vireosnp_dir']}{{num}}/Sample_{{id1}}/{config['donors_vcf']}", zip,  id1=sample_name),
        #expand(f"{config['filt_vcf_dir']}{{num}}/Sample_{{id1}}/{config['filt_vcf']}", zip,  id1=sample_name)
        # check_log_version()




# TO check if the last digit of the line in the error log of STARsolo is a number
def check_isnumber(x):
    try:
        int(x)
        return True
   
    except ValueError:
        return False



# This function reads the log file created per attempt to change the parameter "limitsjdbInsertNsj" in STARsolo   
# We can use this to similarly change other params in the log
def get_limitsjdbval(wildcards, resources, config):
    # This is to check the log file produced after each attempt for the error value
    log_list = glob2.glob("{}{}_STARsolo_log.txt*".format(config['bams_dir'], fold_struct, num=wildcards.num, id1=wildcards.id1))
    for log_file in log_list: 
        with open(log_file) as fin:
            for line in fin:
                if line.startswith("SOLUTION") and "limitSjdbInsertNsj" in line and check_isnumber(line.split()[-1]):
                    return line.split()[-1]
   # This is to check the parameters file (if there was a previous successful run)
    else:
        if os.path.isfile("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)):
            with open("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)) as fin:
                for line in fin:
                    return re.search("--limitSjdbInsertNsj ([0-9]+) ", line).group(1)

        return 1000000



def get_limitsjcollapsed(wildcards, resources, config):
    # This is to check the log file produced after each attempt for the error value
    log_list = glob2.glob("{}{}_STARsolo_log.txt*".format(config['bams_dir'], fold_struct, num=wildcards.num, id1=wildcards.id1))
    for log_file in log_list:
        with open(log_file) as fin:
            for line in fin:
                if line.startswith("SOLUTION") and "limitSjdbInsertNsj" in line and check_isnumber(line.split()[-1]):
                    return line.split()[-1]
           
                elif line.startswith("Solution") and "limitOutSJcollapsed" in line:
                    return 1000000*(1+resources.attempt)

   # This is to check the parameters file (if there was a previous successful run)
    else:
        if os.path.isfile("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)):
            with open("{}Sample_{id1}-cDNA.txt".format(config['star_params_dir'], id1=wildcards.id1)) as fin:
                for line in fin:
                    return int(re.search("--limitOutSJcollapsed ([0-9]+) ", line).group(1))*(1+resources.attempt)

        return 1000000



def get_bams(wildcards, config):
    bam_files = [f for f in glob2.glob("{}{a}/Sample_{b}/{b}_*{c}".format(config['bams_dir'], a=wildcards.num, b=wildcards.id1, c=config['bam']))]
    return bam_files



def get_fastqs(wildcards, config):
    fastq_files = sorted([f for f in glob2.glob("{}{a}/Sample_{b}/{b}_*.R?.fastq.gz".format(config['cDNA_fastqs_dir'], a=wildcards.num, b=wildcards.id1.replace('-cDNA', '-*HTO*')))])
    return fastq_files



def get_STARsolo_mat(wildcards):
    ss_mats = [f for f in glob2.glob("{}{a}/Sample_{b}/{b}_*_Solo.out/GeneFull/filtered_Lun/matrix.mtx.gz".format(config['bams_dir'], a=wildcards.num, b=wildcards.id1))]
    return ss_mats



def get_hto_demux(wildcards):
    hto_d_f = [f for f in glob2.glob("/sc/arion/projects/psychAD/Single_cell_data/preprocess/{a}/HTO_demultiplexed/Sample_{b}_HTOdemux.csv".format(a=wildcards.num, b=wildcards.id1.replace('-cDNA', '-*HTO*')))]
    return hto_d_f



def get_multiseq(wildcards):
    multi_f = [f for f in glob2.glob("/sc/arion/projects/psychAD/Single_cell_data/preprocess/{a}/HTO_demultiplexed/Sample_{b}_MULTIseq.csv".format(a=wildcards.num, b=wildcards.id1.replace('-cDNA', '-*HTO*')))]
    return multi_f



def get_filt_barcodes(wildcards, config):
    barc_f = [f for f in glob2.glob("{}{}{a}_{b}.txt".format(config['filt_barcodes_dir'], config['filt_barcodes'], a=wildcards.num.replace('round_num', ''), b=wildcards.id1))]
    return barc_f



def get_filt_barcodes_cr(wildcards):
    barc_f = [f for f in glob2.glob("/sc/arion/projects/psychAD/pnm/inp_for_vireo/filt_bc_vireo_{a}_{b}_cr.txt".format(a=wildcards.num.replace('round_num', ''), b=wildcards.id1))]
    return barc_f



def calc_donors(wildcards, config):
    temp_df = pd.read_csv(config['wet_lab_info'])
    n_donors = temp_df.loc[temp_df["cDNA_ID"] == wildcards.id1, "SubID"].shape[0]
    return n_donors



def get_donors(wildcards, config):
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
        R1=f"{config['cDNA_fastqs_dir']}{fold_struct}{config['R1_suffix']}",
        R2=f"{config['cDNA_fastqs_dir']}{fold_struct}{config['R2_suffix']}"

    priority: 10

    params:
        gtf=config['gtf_file'],
        genome_dir=config['genome_dir'],
        #output_prefix=lambda wildcards, input: input[0][:-29],
        overhang=config['sjdboverhang'],
        limitsjdbval=get_limitsjdbval,
        chemistry=config['soloType'], # For STARsolo
        whitelist=config['whitelist'], # V3 whitelist
        UMI_length=config['umi_len'], # V3 
        SAM_attr=config['SAM_attr'],
        features=config['features'],
        save_params=f"{config['star_params_dir']}Sample_{{id1}}-cDNA.txt",
        star_def_log_out=f"{config['bams_dir']}{fold_struct}_Log.out",
        limitsjcollap=get_limitsjcollapsed,
        solo_cell_filter=config['solo_cell_filter'],
        out_pref=lambda wildcards, output: output.genefull_lun_matrix[:-13],
        struct_fold=fold_struct
        #cell_filtering not present in 2.7.5b
        #cell_filtering="EmptyDrops_CR" # Cell Filtering matching CellRanger 3.0.0 from Lun et. al. 2019
        # Few more parameters to match CellRanger >=4.0.0

    output:
        f"{config['bams_dir']}{fold_struct}{config['bam']}",
        f"{config['bams_dir']}{fold_struct}{config['STAR_log_final']}",
        f"{config['bams_dir']}{fold_struct}{config['genefull_features']}",
        f"{config['bams_dir']}{fold_struct}{config['genefull_summary']}",
        f"{config['bams_dir']}{fold_struct}{config['gene_features']}",
        f"{config['bams_dir']}{fold_struct}{config['gene_summary']}",
        f"{config['bams_dir']}{fold_struct}{config['barcodes_stats']}",
        f"{config['bams_dir']}{fold_struct}{config['genefull_lun_matrix']}",
        f"{config['bams_dir']}{fold_struct}{config['genefull_lun_features']}",
        f"{config['bams_dir']}{fold_struct}{config['genefull_lun_barcodes']}"


    log:
        f"{config['bams_dir']}{fold_struct}_STARsolo_log.txt"
   
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
        echo "{params.limitsjdbval}, {resources.attempt}"
        if [ ! -d {config[star_params_dir]} ]; then mkdir {config[star_params_dir]}; fi
        STAR --genomeDir {params.genome_dir} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.overhang} --limitSjdbInsertNsj {params.limitsjdbval} --twopassMode Basic --readFilesCommand zcat --readFilesIn {input.R2} {input.R1} --soloType {params.chemistry} --soloUMIlen {params.UMI_length} --soloCBwhitelist {params.whitelist} --soloFeatures {params.features} --soloCellFilter {params.solo_cell_filter} --outSAMattributes {params.SAM_attr} --limitOutSJcollapsed {params.limitsjcollap} --outSAMtype BAM SortedByCoordinate --runThreadN 12 --outFileNamePrefix {config['bams_dir']}{params.struct_fold}_ &> {log}_{resources.attempt}
        gzip {config['bams_dir']}{params.struct_fold}{params.out_pref}*
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
        f"{config['bams_dir']}{fold_struct}{config['bam']}"

    priority: 10

    output:
        f"{config['bams_dir']}{fold_struct}{config['bai']}"

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
        bams=f"{config['bams_dir']}{fold_struct}{config['bam']}"

    priority: 9

    output:
        f"{config['bams_dir']}{fold_struct}{config['gc_bias_metrics']}",
        f"{config['bams_dir']}{fold_struct}{config['gc_summary_metrics']}"

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
        bams=f"{config['bams_dir']}{fold_struct}{config['bam']}"

    priority: 9

    output:
        f"{config['bams_dir']}{fold_struct}{config['rnaseq_metrics']}"

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
        R1="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=fold_struct.replace('-cDNA', '-HTO'), suff=config['R1_suffix']),
        R2="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=fold_struct.replace('-cDNA', '-HTO'), suff=config['R2_suffix'])

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['feature_barcodes']}"

    resources:
        mem_mb=100, #allocate_mem_KBP,
        time_min=10

    shell: 
        """
        python3 helper_py_scripts/create_Feat_Barc.py
        """

         


rule create_mismatch_fasta:
    input:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['feature_barcodes']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['features_mismatch_fa']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['features_mismatch_t2g']}"

    resources:
        mem_mb=100, #allocate_mem_KBP,
        time_min=10


    shell: 
        """
        python3 {config[featuremap_script]} {input} --t2g {output[1]} --fa {output[0]}
        """


rule build_kallisto_index:
    input:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['features_mismatch_fa']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['features_mismatch_idx']}"

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
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['features_mismatch_idx']}",
        R1="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=fold_struct.replace('-cDNA', '-HTO'), suff=config['R1_suffix']),
        R2="{parent_dir}{fs}{suff}".format(parent_dir=config['HTO_fastqs_dir'], fs=fold_struct.replace('-cDNA', '-HTO'), suff=config['R2_suffix'])

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_file']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['ec_matrix']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['tx']}"

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
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_file']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_file_corrected']}"

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
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_file_corrected']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_file_sorted']}"

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
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_file_sorted']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['ec_matrix']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['tx']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['features_mismatch_t2g']}"

    priority: 10

    output:
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_count_dir']}{config['bus_count_mtx']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_count_dir']}{config['bus_count_genes']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_count_dir']}{config['bus_count_barcodes']}"

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
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_count_dir']}{config['bus_count_mtx']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_count_dir']}{config['bus_count_genes']}",
        f"{config['kallisto_bustools_dir']}{fold_struct_kb}{config['bus_count_dir']}{config['bus_count_barcodes']}"

    priority: 10

    output:
        f"{config['h5ad_bustools_dir']}{fold_struct_demux}.h5ad"

    resources:
        mem_mb=3000, #allocate_mem_KBP,
        time_min=10

    shell:
        """
        python3 helper_py_scripts/create_h5ad_from_bustools.py {input[0]} {input[1]} {input[2]}
        """



rule run_calico_solo:
    input:
        f"{config['h5ad_bustools_dir']}{fold_struct_demux}.h5ad",
        starsolo_out=f"{config['bams_dir']}{fold_struct}{config['genefull_lun_matrix']}" #get_STARsolo_mat

    priority: 8
   
    output:
        f"{config['calico_solo_dir']}{fold_struct_demux}{config['calico_solo_h5ad']}"
  
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
        f"{config['calico_solo_dir']}{fold_struct_demux}{config['calico_solo_h5ad']}",
        starsolo_out=f"{config['bams_dir']}{fold_struct}{config['genefull_lun_matrix']}"

    priority: 8

    output:
        f"{config['final_count_matrix_dir']}{fold_struct_demux}{config['final_count_matrix_h5ad']}",
        f"{config['demultiplex_info_dir']}{fold_struct_demux}{config['demultiplex_info']}"

    params:
        mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        min_genes=config['min_genes_per_cell'], # Min #genes per cell
        min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        samples_info=config['wet_lab_info'], # File containing multiplexing info of each set
        genes_info=config['gene_info_file'] # File containing gene names and gene ids for annotations

    
    resources:
        mem_mb=allocate_mem_DXP,
        time_min=30

    shell:
        """
        python3 helper_py_scripts/demul_samples.py {input[0]} {input[1]} {output[0]} {output[1]} -m {params.mito} -g {params.min_genes} -c {params.min_cells}
        sleep 100
        """
        
       


# UMI tag is turned on. Therefore, PCR duplicates are included
rule cellSNP:
    input:
        bc=get_filt_barcodes,
        bams=f"{config['bams_dir']}{fold_struct}{config['bam']}"

    output:
        f"{config['cellsnp_dir']}{fold_struct_kb}{config['cellsnp_cells']}",
        f"{config['cellsnp_dir']}{fold_struct_kb}{config['cellsnp_base']}"

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
        f"{config['cellsnp_dir']}{fold_struct_kb}{config['cellsnp_cells']}",
        config['genotyped_vcf'],
        config['meta_data_geno_samp']
  
    output:
        f"{config['vireosnp_dir']}{fold_struct_kb}{config['donors_vcf']}",
        f"{config['filt_vcf_dir']}{fold_struct_kb}{config['filt_vcf']}"

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
        bam=f"{config['bams_dir']}{fold_struct}{config['bam']}",
        barcode_list="" # Barcodes vs sample name txt file, produced after producing final_count_matrices

    params:
        output_pref=f"{config['split_bams_dir']}{fold_struct}_"
    output:
        f"{config['split_bams_proxy_dir']}{{num}}_Sample_{{id1}}.txt" # Proxy to the output

    shell:
        """
        helper_py_scripts/split_bam_indiv_barc_samtools.sh {input.bam} {input.barcode_list} {params.output_pref} > {output}
        sleep 100
        """



# Rule to run update.py
rule update_stats:
    input:
        STAR_log=expand(f"{config['bams_dir']}{fold_struct}{config['STAR_log_final']}", zip,  id1=sample_name),
        PICARD_GC=expand(f"{config['bams_dir']}{fold_struct}{config['gc_summary_metrics']}", zip,  id1=sample_name),
        PICARD_RNAseq=expand(f"{config['bams_dir']}{fold_struct}{config['rnaseq_metrics']}", zip,  id1=sample_name),
        SS_G_Feat=expand(f"{config['bams_dir']}{fold_struct}{config['gene_features']}", zip,  id1=sample_name),
        SS_GF_Feat=expand(f"{config['bams_dir']}{fold_struct}{config['genefull_features']}", zip,  id1=sample_name),
        SS_G_Summ=expand(f"{config['bams_dir']}{fold_struct}{config['gene_summary']}", zip,  id1=sample_name),
        SS_GF_Summ=expand(f"{config['bams_dir']}{fold_struct}{config['genefull_summary']}", zip,  id1=sample_name),
        SS_Barcodes=expand(f"{config['bams_dir']}{fold_struct}{config['barcodes_stats']}", zip,  id1=sample_name),
        Demultiplex_info=expand(f"{config['demultiplex_info_dir']}{fold_struct_demux}{config['demultiplex_info']}", zip,  id1=sample_name)

    params:
        map_file=config['meta_data']

    output:
        config['log_all_stats']

    threads: 1

    resources:
        mem_mb=500,
        time_min=30

    script:
        "helper_py_scripts/update_logs.py"


#rule 
