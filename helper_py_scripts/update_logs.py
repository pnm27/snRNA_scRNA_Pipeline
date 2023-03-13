#!/usr/bin/env python3

import pandas as pd
import os
import itertools
import numpy as np
import glob2
import re
from time import sleep
import errno, argparse


# This function returns the filename if it exists
# otherwise an empty string
def get_filename(loc_dir, file_struct, fn, suffix):
    if file_struct.endswith('/'):
        try:
            return glob2.glob(os.path.join(loc_dir, file_struct, f"{fn}*{suffix}"))[0]
        except:
            return ""
    elif file_struct == "":
        try:
            return glob2.glob(os.path.join(loc_dir, f"{fn}*{suffix}"))[0]
        except:
            return ""
    else:
        if glob2.glob(os.path.join(loc_dir, f"{file_struct}*{suffix}")):
            return glob2.glob(os.path.join(loc_dir, f"{file_struct}*{suffix}"))[0]
        elif glob2.glob(os.path.join(loc_dir, f"{file_struct}*{fn}*{suffix}")):
            return glob2.glob(os.path.join(loc_dir, f"{file_struct}*{fn}*{suffix}"))[0]
        else:
            return ""



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


# Main function to write rows of samples into a df
def write_logs(big_df, mapper, all_files_dict, no_progs, **kwargs):
   
    new_row = []
    # Extra Annotations through kwargs
    # Reproduce the sequence of how the new_columns were set
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
        add_value=""
        if prog != "LAB" and sub_prog not in no_progs:
            if sub_prog == "REG":
                temp_df = pd.read_csv(all_files_dict["STAR_final"], names=["cols", "vals"], delimiter=r"|", skiprows=[7, 22, 27, 34])
                temp_df["vals"] = temp_df.vals.str.strip()
                temp_df["cols"] = temp_df.cols.str.strip()
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[mapper["curr_val"] == val, "val_in_log"].values[0], "vals"].values[0]
                except:
                    add_value = ""
                new_row.append(add_value.replace(" ","/"))

            elif sub_prog == "GC":
                temp_df = pd.read_csv(all_files_dict["PICARD_GC"], sep='\t', skiprows=6)
                try:
                    add_value = temp_df.loc[0, mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GC"), "val_in_log"].values[0]]
                except:
                    add_value = ""
                new_row.append(add_value)
               
            elif sub_prog == "RNASEQMETRIC":
                temp_df = pd.read_csv(all_files_dict["PICARD_RNASeq"], sep='\t', nrows=1, skiprows=6)
                try:
                    add_value = temp_df.loc[0, mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "RNASEQMETRIC"), "val_in_log"].values[0]]
                except:
                    add_value = ""
                new_row.append(add_value)

            elif sub_prog == "GENE_FEATURE":
                temp_df = get_df(all_files_dict["Gene_Features"]) 
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENE_FEATURE"), "val_in_log"].values[0], "vals"].values[0]
                except:
                    add_value = ""
                new_row.append(add_value)

            elif sub_prog == "GENE_SUMM":
                temp_df = pd.read_csv(all_files_dict["Gene_Summary"], names=['cols', 'vals'])
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENE_SUMM"), "val_in_log"].values[0], "vals"].values[0]
                except:
                    add_value = ""
                new_row.append(add_value)

            elif sub_prog == "GENEFULL_FEATURE":
                temp_df = get_df(all_files_dict["GeneFull_Features"])
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENEFULL_FEATURE"), "val_in_log"].values[0], "vals"].values[0]
                except:
                    add_value = ""
                new_row.append(add_value)

            elif sub_prog == "GENEFULL_SUMM":
                temp_df = pd.read_csv(all_files_dict["GeneFull_Summary"], names=['cols', 'vals'])
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "GENEFULL_SUMM"), "val_in_log"].values[0], "vals"].values[0]
                except:
                    add_value = ""
                new_row.append(add_value)

            elif sub_prog == "BARCODE_STATS":
                temp_df = get_df(all_files_dict["Barcodes_stats"])
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "BARCODE_STATS"), "val_in_log"].values[0], "vals"].values[0]
                except:
                    add_value = ""
                new_row.append(add_value)
           
            elif sub_prog == "DEMUX":
                temp_df = pd.read_csv(all_files_dict["Demultiplex_stats"], names=['cols', 'vals'], skiprows=1, sep='\t')
                add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "DEMUX"), "val_in_log"].values[0], "vals"].values[0]
                if val == "N_CELLS_AFTER_DEMUX_CS" and add_value.endswith(','):
                    new_row.append(add_value[:-1])

                # Compatibility with older-style of producing demux_info file
                elif val == "N_CELLS_AFTER_DEMUX_CS" and 'Name:' in add_value:
                   add_value = add_value[:add_value.find('Name:')]
                   add_value= re.sub('[^\S\r\n]+', ':', add_value)
                   add_value = add_value[:-1]
                   add_value = re.sub('\n', ',', add_value)
                   new_row.append(add_value)

                else:
                     new_row.append(add_value)         

            else:
                raise ValueError(f'This extra column exists in the output file-All_logs.csv: {prog}, {sub_prog}, {val}')

        elif prog != "LAB" and sub_prog in no_progs:
            new_row.append(add_value)

        else:
            continue


    return new_row


# Extra columns (annotations) to add (This sequence should be maintained everywhere)
new_cols_to_add = [['ROUND', 'LAB', 'BATCH'], ['SAMPLE', 'LAB', 'SAMPLE'], ['SET', 'LAB', 'BATCH'], ['PREPARER', 'LAB', 'BATCH'], ['REP', 'LAB', 'BATCH']]

# Function to conditionally run this script through Snakemake if the current file has fewer columns that the last version of this script
def get_latest_extra_columns():
    global new_cols_to_add
    return len(new_cols_to_add)

def get_argument_parser():
    """Generate and return argument parser."""

    parser = argparse.ArgumentParser(description="Compile all the files to combine all stats. NOTE: For all optional files (including all folder structures), if parameter is present but no value is \
        provided then values will be used as described by the defaults in respective help message. Folder structures for STARsolo output and PICARD is assumed to be the same, by default. \
        For optional files, the only required value is that of the parent folder.")

    parser.add_argument('samples', nargs='+', help="List of samples. This will be the prefix(es) for all the name(s) of the output file(s)")

    
    # Optional parameters
    parser.add_argument('-m', '--map_file', help="Mapping file that contains info on the headers in the output. DEFAULT: <current working directory>/Final_out_MAP_2.tsv", 
        default=os.path.join(os.getcwd()+"Final_out_MAP_2.tsv"))
    parser.add_argument('-o', '--output_file', help="output file. DEFAULT: <current working directory>/All_logs.tsv", default=os.path.join(os.getcwd()+"All_logs.tsv"))
    parser.add_argument('-b', '--bam_dir', help="Directory containing bam file(s). DEFAULT: current working directory", default=os.getcwd())
    parser.add_argument('-p', '--picard_dir', help="Directoy containing PICARD outputs. DEFAULT: current working directory", default=os.getcwd())
    parser.add_argument('-d', '--demul_dir', help="Directory containing demultiplexing stats. DEFAULT: current working directory", default=os.getcwd())
    parser.add_argument('--bam_struct', help="Regex to identify bam file(s) for the give sample(s). NOTE: In the regex, <sample> denotes where to insert the sample name(s) \
        provided to this script. DEFAULT: \"<current_working_dir>/<sample>/\"", default=os.path.join(os.getcwd(), "<sample>/"))
    parser.add_argument('--pc_struct', help="Regex to identify picard file(s) for the give sample(s). NOTE: In the regex, <sample> denotes where to insert the sample name(s) \
        provided to this script. DEFAULT: \"<current_working_dir>/<sample>/\"", default=os.path.join(os.getcwd(), "<sample>/"))
    parser.add_argument('--dem_struct', help="Regex to identify demultiplex info containing file(s) for the give sample(s). NOTE: In the regex, <sample> denotes where to insert the sample name(s) \
        provided to this script. DEFAULT: \"<current_working_dir>/<sample>\"", default=os.path.join(os.getcwd(), "<sample>"))
    parser.add_argument('--ss_l', nargs='?', help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_Log.final.out\"", const="_Log.final.out", default=None)
    parser.add_argument('--pc_gc', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_summary_metrics.txt\"", const="_summary_metrics.txt", default=None)
    parser.add_argument('--pc_rs', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_rnaseq_metrics.txt\"", const="_rnaseq_metrics.txt", default=None)
    parser.add_argument('--ss_g_f', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_Solo.out/Gene/Features.stats\"", const="_Solo.out/Gene/Features.stats", default=None)
    parser.add_argument('--ss_gf_f', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_Solo.out/GeneFull/Features.stats\"", const="_Solo.out/GeneFull/Features.stats", default=None)
    parser.add_argument('--ss_g_s', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_Solo.out/Gene/Summary.csv\"", const="_Solo.out/Gene/Summary.csv", default=None)
    parser.add_argument('--ss_gf_s', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_Solo.out/GeneFull/Summary.csv\"", const="_Solo.out/GeneFull/Summary.csv", default=None)
    parser.add_argument('--ss_bc', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_Solo.out/Barcodes.stats\"", const="_Solo.out/Barcodes.stats", default=None)
    parser.add_argument('--dem_info', nargs='?',  help="Suffix for the output (if not the same as the default one). Absence of this parameter is treated as not intended in the compilation. \
        DEFAULT: \"_STARsolo_info.tsv\"", const="_STARsolo_info.tsv", default=None)

    return parser


def main():
    """Main entry point"""

    # Parse arguments
    parser = get_argument_parser()
    args = parser.parse_args()
    
    # Parse directory values
    bam_dir = args.bam_dir
    pic_dir = args.picard_dir
    dem_dir = args.demul_dir



    # Validate the optional parameters, if present
    opt_file_params = {'ss_l': ['.out', 'REG'], 'pc_gc': ['.txt', 'GC'], 'pc_rs':['.txt', 'RNASEQMETRIC'], 'ss_g_f': ['.stats', 'GENE_FEATURE'], 'ss_gf_f': ['.stats', 'GENEFULL_FEATURE'],
                       'ss_g_s': ['.csv', 'GENE_SUMM'], 'ss_gf_s': ['.csv', 'GENEFULL_SUMM'], 'ss_bc': ['.stats', 'BARCODE_STATS'], 'dem_info': ['.tsv', 'DEMUX'], 'output_file': ['.tsv', None],
                       'map_file': ['.tsv', None]}

    # List of programs from which no stats need be recorded
    exclude_progs=[]

    # Extension test for files and dir exists for directories
    for k, v in vars(args).items():
        # If these parameters are present, they should have appropriate extensions
        if k in opt_file_params and v != None and not v.endswith(opt_file_params[k][0]):
            raise ValueError(f"The file extension in {v} for the parameter {k} is unexpected!")
        elif k.endswith('dir') and ( v == None or not os.path.isdir(v)):
            raise ValueError(f"The directory {v} provided for the parameter {k} doesn't exist!")
        elif k in opt_file_params and v == None and not k.endswith('dir') and not k.endswith('file'):
            exclude_progs.append(opt_file_params[k][1])
        else:
            continue


    out=args.output_file
    map_names = pd.read_csv(args.map_file, delimiter="\t", names=["val_in_log", "curr_val", "prog", "sub_prog", "desc"])
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



    # Process each sample -----------------------------------------------------------------------------------------------------------------------------------------------
    # List containing per sample values as lists (list of lists)
    row_list = []
    for sample in args.samples:


        # create per sample copy of exclude_prog list
        samp_excl_progs = exclude_progs.copy()

        # Parse file structures for bam files, picard files and demultiplex info files
        bam_st = args.bam_struct.replace("<sample>", sample)
        pc_st = args.pc_struct.replace("<sample>", sample)
        dem_st = args.dem_struct.replace("<sample>", sample)
        
        # print("Entered loop")
        # Get full filenames if user requires them to be tabulated
        ss_log_final = get_filename(bam_dir, bam_st, sample, args.ss_l) if args.ss_l != None else ""
        ss_gene_summary = get_filename(bam_dir, bam_st, sample, args.ss_g_s) if args.ss_g_s != None else ""
        ss_genefull_summary = get_filename(bam_dir, bam_st, sample, args.ss_gf_s) if args.ss_gf_s != None else ""
        ss_gene_features = get_filename(bam_dir, bam_st, sample, args.ss_g_f) if args.ss_g_f != None else ""
        ss_genefull_features = get_filename(bam_dir, bam_st, sample, args.ss_gf_f) if args.ss_gf_f != None else ""
        ss_bc_stats = get_filename(bam_dir, bam_st, sample, args.ss_bc) if args.ss_bc != None else ""
        pc_gc_file = get_filename(pic_dir, pc_st, sample, args.pc_gc) if args.pc_gc != None else ""
        pc_rs_file = get_filename(pic_dir, pc_st, sample, args.pc_rs) if args.pc_rs != None else ""
        dem_file = get_filename(dem_dir, dem_st, sample, args.dem_info) if args.dem_info != None else ""

        # Example for Adding additional info per sample
        sample_name=sample
        # r_num = int(re.search('/round([0-9]+)/', files_dict["STAR_final"]).group(1))
        preparer = sample.split('-')[2][0]
        replicate = sample.split('-')[2][1]
        set_val = sample[:-6]

        
        # Test if at least one of the input files exists
        test_f_exists=[ss_log_final == "", ss_gene_summary == "", ss_genefull_summary == "", ss_gene_features == "", ss_genefull_features == "", ss_bc_stats == "",
            pc_gc_file == "", pc_rs_file == "", dem_file == ""]
            
        if not any(test_f_exists):
            raise ValueError(f"All files for the sample {sample} are empty or not to be found! Please check the directories and usage of this script for more info.")


        files_dict = {"STAR_final": ss_log_final, "PICARD_GC": pc_gc_file, "PICARD_RNASeq": pc_rs_file, "Gene_Features": ss_gene_features, "GeneFull_Features": ss_genefull_features, 
                      "Gene_Summary": ss_gene_summary, "GeneFull_Summary": ss_genefull_summary, "Barcodes_stats": ss_bc_stats, "Demultiplex_stats": dem_file}


        # Check for each sample in the list if it has all required files otherwise mark as "" for the respective sample (i.e. update samp_excl_progs list)
        per_samp_check = {'REG': ss_log_final, 'GC': pc_gc_file, 'RNASEQMETRIC': pc_rs_file, 'GENE_FEATURE': ss_gene_features, 'GENEFULL_FEATURE': ss_genefull_features,
                        'GENE_SUMM': ss_gene_summary, 'GENEFULL_SUMM': ss_genefull_summary, 'BARCODE_STATS': ss_bc_stats, 'DEMUX': dem_file}
        for k, v in per_samp_check.items():
            if k not in samp_excl_progs and v == "":
                samp_excl_progs.append(k)


           
        if not(combo_log['LAB']['SAMPLE']['SAMPLE'].str.contains(sample).any()) :
            # Add a kwargs style input for extra annotations
            row_list.append(write_logs(combo_log, map_names, files_dict, samp_excl_progs, sample=sample_name, prep=preparer, rep=replicate, set_num=set_val))

        print(f"Finished adding {sample_name} to the file")


    temp_df = pd.DataFrame(row_list, columns=pd.MultiIndex.from_frame(cl, names=["prog", "sub_prog", "curr_val"]))
    combo_log = pd.concat([combo_log, temp_df], ignore_index=True)


    # If the demux file is not used for compilation then skip these steps
    if args.dem_info != None:
        try:
            # combo_log[("STARsolo", "DEMUX", "DOUBLET_PCT")] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_DOUBLET_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
            combo_log[("STARsolo", "DEMUX", "DOUBLET_PCT")] = combo_log[("STARsolo", "DEMUX", "N_DOUBLET_CELLS_CS")].astype(int)/combo_log[("STARsolo", "DEMUX", "N_CELLS_START")].astype(int)
        except:
            print("Can't calculate Doublet ratio! Check output file for more info!")

        try:
            # combo_log[("STARsolo", "DEMUX", "NEGATIVE_PCT")] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_NEGATIVE_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
            combo_log[("STARsolo", "DEMUX", "NEGATIVE_PCT")] = combo_log[("STARsolo", "DEMUX", "N_NEGATIVE_CELLS_CS")].astype(int)/combo_log[("STARsolo", "DEMUX", "N_CELLS_START")].astype(int)
        except:
            print("Can't calculate Negative ratio! Check output file for more info!")

        try:
            combo_log[("STARsolo", "DEMUX", "N_DEMUXED_CELLS")] = (combo_log[("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT")].astype(int)-combo_log[("STARsolo", "DEMUX", "N_DOUBLET_CELLS_CS")].astype(int).values
                        -combo_log[("STARsolo", "DEMUX", "N_NEGATIVE_CELLS_CS")].astype(int).values)

        except:
            print("Can't calculate percentage of cells retained for demultiplexing")

        try:
            combo_log[("STARsolo", "DEMUX", "CELL_RENTENTION")] = combo_log[("STARsolo", "DEMUX", "N_DEMUXED_CELLS")]/combo_log[("STARsolo", "DEMUX", "N_CELLS_START")].astype(int).values
        except:
            print("Can't calculate percentage of cells retained for demultiplexing")


    combo_log.to_csv(out, sep = "\t", index=False)


    sleep(30)
    

# Run this only when executed through Snakemake
if __name__ == "__main__":
    main()
