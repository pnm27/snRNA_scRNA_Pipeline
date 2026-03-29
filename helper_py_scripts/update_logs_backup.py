#!/usr/bin/env python3

import pandas as pd
import os, json
import numpy as np
import glob2
from time import sleep
import argparse, warnings # errno
from demultiplex_helper_funcs import process_columns, auto_read
from jsonschema import validate, Draft202012Validator


# This function returns the filename if it exists
# otherwise an empty string
def get_filename(loc_dir, file_struct, fn, suffix):

    if suffix is None:
        return ''
    
    loc_dir = [loc_dir] if ',' not in loc_dir else loc_dir.split(',')
    ret_path = []
    for l in loc_dir:
        if file_struct.endswith('/'):
            try:
                ret_path.append(glob2.glob(os.path.join(l, file_struct, f"{fn}*{suffix}"))[0])
            except:
                ret_path.append("")
        elif file_struct == "":
            try:
                ret_path.append(glob2.glob(os.path.join(l, f"{fn}*{suffix}"))[0])
            except:
                ret_path.append("")
        else:
            if glob2.glob(os.path.join(l, f"{file_struct}*{suffix}")):
                ret_path.append(glob2.glob(os.path.join(l, f"{file_struct}*{suffix}"))[0])
            elif glob2.glob(os.path.join(l, f"{file_struct}*{fn}*{suffix}")):
                ret_path.append(glob2.glob(os.path.join(l, f"{file_struct}*{fn}*{suffix}"))[0])
            else:
                ret_path.append("")

    ret_path = ','.join(ret_path)
    if ret_path == ',':
        return ''
    else:
        return ret_path



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
def write_logs(output_df_cols, mapper, all_files_dict, no_progs, more_data): # had **kwargs
   
    new_row = []
    # Extra Annotations through kwargs
    # Reproduce the sequence of how the new_columns were set
    # OLD STYLE
    # new_row.append(kwargs["sample"])
    # new_row.append(kwargs["set_num"])
    # new_row.append(kwargs["prep"])
    # new_row.append(kwargs["rep"])
    new_row.extend([ f.column_value for f in more_data ])
    # NEW STYLE
    

    # Store Barcode.stats and Feature.stats files as DF in a list for easy access
    # Seq is Feature.stats for Gene, Feature.stats for GeneFull, Barcode.stats
    # stats_file = [x for k, x in all_files_dict.items() if '.stats' in x ]
    # list_df = [get_df(x) for x in stats_file ]
    
    
    # Add values to a list in the same sequence as the final output file/dataframe
    for prog, sub_prog, val in output_df_cols:
        add_value=""
        if prog != "LAB" and (sub_prog not in no_progs and
                        not any(list(map(lambda x: sub_prog.startswith(x), no_progs)))
                        ):
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
           
            elif sub_prog.startswith("DEMUX"):
                temp_df = pd.read_csv(all_files_dict["Demultiplex_stats"], names=['cols', 'vals'], skiprows=1, sep='\t')
                try:
                    add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == sub_prog), "val_in_log"].values[0], "vals"].values[0]
                    if isinstance(add_value, str) and add_value.endswith(','):
                        new_row.append(add_value[:-1])
                    elif isinstance(add_value, float):
                        new_row.append(add_value if float(add_value).is_integer else round(add_value, 3))
                    else:
                        new_row.append(add_value)
                except:
                    # print("FIND VALUE: ")
                    # print(f"current_val: {val}, sub_prog: {sub_prog}")
                    add_value = ""
                    new_row.append(add_value)
                    
                # DEPRACATED WITH ADVENT OF NEW SUB_PROG 'DEMUX_CS' AND 'DEMUX_VS'
                # if val == "N_CELLS_AFTER_DEMUX_CS" and add_value.endswith(','):
                #     new_row.append(add_value[:-1])

                # # Compatibility with older-style of producing demux_info file
                # elif val == "N_CELLS_AFTER_DEMUX_CS" and 'Name:' in add_value:
                #    add_value = add_value[:add_value.find('Name:')]
                #    add_value= re.sub('[^\S\r\n]+', ':', add_value)
                #    add_value = add_value[:-1]
                #    add_value = re.sub('\n', ',', add_value)
                #    new_row.append(add_value)

                # else:
                #     new_row.append(add_value)


            # elif sub_prog == "DEMUX_CS":
            #     temp_df = pd.read_csv(all_files_dict["Demultiplex_stats"], names=['cols', 'vals'], skiprows=1, sep='\t')
            #     add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "DEMUX_CS"), "val_in_log"].values[0], "vals"].values[0]
            #     if add_value.endswith(','):
            #         new_row.append(add_value[:-1])
            #     else:
            #          new_row.append(add_value)        

            # elif sub_prog == "DEMUX_VS":
            #     temp_df = pd.read_csv(all_files_dict["Demultiplex_stats"], names=['cols', 'vals'], skiprows=1, sep='\t')
            #     add_value = temp_df.loc[temp_df["cols"] == mapper.loc[(mapper["curr_val"] == val) & (mapper["sub_prog"] == "DEMUX_VS"), "val_in_log"].values[0], "vals"].values[0]
            #     if add_value.endswith(','):
            #         new_row.append(add_value[:-1])
            #     else:
            #          new_row.append(add_value)          

            else:
                raise ValueError(f'This extra column exists in the output file-All_logs.csv: {prog}, {sub_prog}, {val}')

        elif prog != "LAB" and not (sub_prog not in no_progs and
                        not any(list(map(lambda x: sub_prog.startswith(x), no_progs)))
                        ):
            new_row.append("")

        else:
            continue


    return new_row

# OLD-STYLE
# Extra columns (annotations) to add (This sequence should be maintained everywhere)
# new_cols_to_add = [
#     ['ROUND', 'LAB', 'BATCH'], 
#     ['SAMPLE', 'LAB', 'SAMPLE'], 
#     ['SET', 'LAB', 'BATCH'], 
#     ['PREPARER', 'LAB', 'BATCH'], 
#     ['REP', 'LAB', 'BATCH']
# ]


# NOT NEEDED
# Function to conditionally run this script through Snakemake if the current file has fewer columns that the last version of this script
# def get_latest_extra_columns():
#     global new_cols_to_add
#     return len(new_cols_to_add)

def get_all_columns(map_df: pd.DataFrame, new_cols_to_add: list[list[str]]):
    cols = map_df.iloc[:, 1:-1]
    if new_cols_to_add:
        cl = pd.DataFrame(new_cols_to_add, columns=list(map_df.columns.values)[1:-1])
        # cl = cl.append(cols, ignore_index=True) # DEPRACATED
        cl = pd.concat([cl, cols], ignore_index=True, axis=0)

        # Change column header order to: curr_val, sub_prog, prog so that the output file looks like:
        #Prog
        #Sub_prog
        #curr_val

        #eg:
        #STARsolo
        #GENE_SUMM
        #N_READS
        return cl.iloc[:, [1, 2, 0]]
    else:
        return cols.iloc[:, [1, 2, 0]]



def get_argument_parser():
    """Generate and return argument parser."""

    parser = argparse.ArgumentParser(description="Compile all the files to combine all stats. NOTE: For all optional files (including all folder structures), if parameter is present but no value is \
        provided then values will be used as described by the defaults in respective help message. Folder structures for STARsolo output and PICARD is assumed to be the same, by default. \
        For optional files, the only required value is that of the parent folder.")

    # OLD STYLE
    # parser.add_argument('samples', nargs='+', help="List of samples. This will be the prefix(es) for all the name(s) of the output file(s)")
    # New preferred argument
    parser.add_argument(
        "--samples", nargs='+',
        help="Preferred: List of samples or File (comma separated and with headers). This will be the prefix(es) for all the name(s) of the output file(s). \
            describing which version of demultiplexing result to use per each pool. Absence of this \
            parameter is treated as each pool fetches data from the same dir"
    )

    # Old deprecated positional argument
    parser.add_argument(
        "samples_deprecated",
        nargs="*",  # makes it optional
        help=argparse.SUPPRESS  # hide from help output
    )

    # Optional parameters
    parser.add_argument('-m', '--map_file', help="Mapping file that contains info on the headers in the output. DEFAULT: <current working directory>/Final_out_MAP_2.tsv", 
        default=os.path.join(os.getcwd()+"Final_out_MAP_2.tsv"))
    parser.add_argument('-o', '--output_file', help="output file. DEFAULT: <current working directory>/All_logs.tsv", default=os.path.join(os.getcwd()+"All_logs.tsv"))
    parser.add_argument('-b', '--bam_dir', help="Directory containing bam file(s). DEFAULT: current working directory", default=os.getcwd())
    parser.add_argument('-p', '--picard_dir', nargs='?', help="Directoy containing PICARD outputs. DEFAULT: current working directory", 
        const=os.getcwd(), default=None)
    parser.add_argument('-d', '--demul_dir', nargs='?', help="Directory (comma-sep) containing demultiplexing stats. DEFAULT: current working directory", 
        const=os.getcwd(), default=None)
    parser.add_argument('--bam_struct', help="Regex to identify bam file(s) for the give sample(s). NOTE: In the regex, <sample> denotes where to insert the sample name(s) \
        provided to this script. DEFAULT: \"<current_working_dir>/<sample>/\"", default=os.path.join(os.getcwd(), "<sample>/"))
    parser.add_argument('--pc_struct', nargs='?', help="Regex to identify picard file(s) for the give sample(s). NOTE: In the regex, <sample> denotes where to insert the sample name(s) \
        provided to this script. DEFAULT: \"<current_working_dir>/<sample>/\"",const=os.path.join(os.getcwd(), "<sample>/"), 
        default=None)
    parser.add_argument('--dem_struct', nargs='?', help="Regex to identify demultiplex info containing file(s) for the give sample(s). NOTE: In the regex, <sample> denotes where to insert the sample name(s) \
        provided to this script. DEFAULT: \"<current_working_dir>/<sample>\"", const=os.path.join(os.getcwd(), "<sample>"),
        default=None)
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
    parser.add_argument('--common_annotations', nargs='?', help="json file that has extra annotations. By default, it will get these options from the file called annotate_h5ad.json in the root folder. \
    They will be validated against the schema present in the schema dir.", const="demul_samples_annotate_h5ad.json", default=None
    )
    parser.add_argument('-w', '--wet_lab_file', help="Path to file that contains either/all of: HTO info for each set, annotations, etc."
    )

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

    # Resolve which one to use
    if args.samples and args.samples_deprecated:
        parser.error("Provide either --samples or positional samples, not both.")

    if args.samples_deprecated:
        warnings.warn(
            "Positional 'samples' is deprecated. Use --samples instead.",
            DeprecationWarning
        )

    samples = args.samples if args.samples else args.samples_deprecated

    # If annotations is needed then provide the wet lab file
    if args.common_annotations is not None and args.wet_lab_file is None:
        warnings.warn(
            "Annotations are expected as json was provided but Wet lab file was NOT provided. \
            No extra annotations will be added!",
            UserWarning
        )

    # If arguments -p (picard_dir) and -d (demultiplex_dir) is not provided
    # then pc_rs and pc_gc; dem_info are implied to be absent
    # i.e. without a picard output dir there can't be outputs for
    # CollectRnaSeqMetrics and CollectGcBiasMetrics similarly for demultiplexing
    dir_depend_args = {
        'picard_dir': ['pc_gc', 'pc_rs'],
        'demul_dir': ['dem_info'],
    }

    # Validate the optional parameters, if present
    opt_file_params = {
                    'ss_l': ['.out', 'REG'], 'pc_gc': ['.txt', 'GC'], 
                    'pc_rs':['.txt', 'RNASEQMETRIC'], 
                    'ss_g_f': ['.stats', 'GENE_FEATURE'], 
                    'ss_gf_f': ['.stats', 'GENEFULL_FEATURE'],
                    'ss_g_s': ['.csv', 'GENE_SUMM'], 
                    'ss_gf_s': ['.csv', 'GENEFULL_SUMM'], 
                    'ss_bc': ['.stats', 'BARCODE_STATS'], 
                    'dem_info': ['.tsv', 'DEMUX'], 
                    'output_file': ['.tsv', None],
                    'map_file': ['.tsv', None]
                    }

    # List of programs from which no stats need be recorded
    exclude_progs=[]

    # Extension test for files and dir exists for directories (only for BAM)
    # If these parameters are present, they should have appropriate extensions
    for k, v in vars(args).items():
        if k in opt_file_params and v != None and not v.endswith(opt_file_params[k][0]):
            raise ValueError(f"The file extension in {v} for the parameter {k} is unexpected!")
        elif k == 'bam_dir' and ( v == None or not os.path.isdir(v)):
            raise ValueError(f"The directory {v} provided for the parameter {k} doesn't exist!")
        elif k.endswith('dir') and ( v == None or not os.path.isdir(v)):
            print(f"No directory was provided for the parameter {k}! Hence, will remove corresponding\n"
                  " outputs.")
            for p in dir_depend_args[k]:
                exclude_progs.append(opt_file_params[p][1])
        elif k in opt_file_params and v == None and not k.endswith('dir') and not k.endswith('file'):
            exclude_progs.append(opt_file_params[k][1])
        else:
            continue


    out=args.output_file
    map_names = pd.read_csv(args.map_file, delimiter="\t", names=["val_in_log", "curr_val", "prog", "sub_prog", "desc"])

    # If file doesn't exist create one else open as pandas dataframe
    try:
        #if os.path.isfile(snakemake.output[0]) :
        combo_log = pd.read_csv(out, sep = "\t", header=[0, 1, 2])

        # DON'T SUPPORT CATCHING OLDER FILES
        # Catch older files that may have lesser columns than expected
        # if combo_log.shape[1] != cl.shape[0]:
        #     combo_log = pd.DataFrame(columns=pd.MultiIndex.from_frame(cl, names=["prog", "sub_prog", "curr_val"]))


    except:
        #with open(snakemake.output[0], 'w+') as fout:
        # combo_log = pd.DataFrame(columns=pd.MultiIndex.from_frame(cl, names=["prog", "sub_prog", "curr_val"]))
        combo_log = []

    if args.common_annotations is not None:
        # Wet Lab file, Filter wet lab file's columns, if needed
        df = auto_read(args.wet_lab_file)

        # Load the user config
        with open(args.common_annotations) as f:
            data = json.load(f)

        # Load your schema (replace with your actual JSON schema file)
        with open("../schema/demul_samples_annotate_h5ad.schema.json") as f:
            schema = json.load(f)
        # Validate
        validator = Draft202012Validator(schema)
        errors = sorted(validator.iter_errors(data), key=lambda e: e.path)

        if errors:
            for err in errors:
                print("Validation error:", err.message)
        else:
            print("Validation successful!")

    # Get all columns for the output df
    extra_cols = [ col.get("columnInLogs", []) for col in data["columns"] ]
    cols_list = get_all_columns(map_names, extra_cols)

    # Process each sample -----------------------------------------------------------------------------------------------------------------------------------------------
    # List containing per sample values as lists (list of lists)
    row_list = []
    for sample in samples:
        skip_sample = False # By default, don't skip any sample

        # create per sample copy of exclude_prog list
        samp_excl_progs = exclude_progs.copy()

        # Parse file structures for bam files, picard files and demultiplex info files
        bam_st = args.bam_struct.replace("<sample>", sample)
        pc_st = args.pc_struct.replace("<sample>", sample) if args.pc_struct is not None else ""
        dem_st = args.dem_struct.replace("<sample>", sample) if args.dem_struct is not None else ""
        
        # print("Entered loop")
        # Get full filenames if user requires them to be tabulated

        ss_log_final = get_filename(bam_dir, bam_st, sample, args.ss_l)
        ss_gene_summary = get_filename(bam_dir, bam_st, sample, args.ss_g_s)
        ss_genefull_summary = get_filename(bam_dir, bam_st, sample, args.ss_gf_s)
        ss_gene_features = get_filename(bam_dir, bam_st, sample, args.ss_g_f)
        ss_genefull_features = get_filename(bam_dir, bam_st, sample, args.ss_gf_f)
        ss_bc_stats = get_filename(bam_dir, bam_st, sample, args.ss_bc)
        pc_gc_file = get_filename(pic_dir, pc_st, sample, args.pc_gc)
        pc_rs_file = get_filename(pic_dir, pc_st, sample, args.pc_rs)
        dem_file = get_filename(dem_dir, dem_st, sample, args.dem_info)

        # sample_name=sample
        # DEPRECATED: Add Annotations
        # preparer = sample.split('-')[2][0] if len(sample.split('-')[2]) > 1 else "NA"
        # replicate = sample.split('-')[2][1] if len(sample.split('-')[2]) > 1 else sample.split('-')[2][0]
        # set_val = '-'.join(sample.split('-')[:2]) + '-' + preparer if preparer != 'NA' else '-'.join(sample.split('-')[:2])
        # NEW STYLE ANNOTIONS
        # PROCESS THE DATA
        processed_data = process_columns(data, sample, df)

        
        # Test if at least one of the input files exists
        test_f_exists = [
            ss_log_final != "", ss_gene_summary != "", 
            ss_genefull_summary != "", ss_gene_features != "", 
            ss_genefull_features != "", ss_bc_stats != "",
            pc_gc_file != "", pc_rs_file != "", dem_file != ""
        ]
            
        if not any(test_f_exists):
            raise ValueError(
                f"All files for the sample {sample} are "
                "empty or not to be found! Please check the directories "
                "and usage of this script for more info."
            )

        files_dict = {
            "STAR_final": ss_log_final, "PICARD_GC": pc_gc_file, 
            "PICARD_RNASeq": pc_rs_file, "Gene_Features": ss_gene_features, 
            "GeneFull_Features": ss_genefull_features, "Gene_Summary": ss_gene_summary, 
            "GeneFull_Summary": ss_genefull_summary, "Barcodes_stats": ss_bc_stats, 
            "Demultiplex_stats": dem_file
        }

        # Check for each sample in the list if it has all required files otherwise mark as "" for the respective sample (i.e. update samp_excl_progs list)
        per_samp_check = {
            'REG': ss_log_final, 'GC': pc_gc_file, 
            'RNASEQMETRIC': pc_rs_file, 'GENE_FEATURE': ss_gene_features, 
            'GENEFULL_FEATURE': ss_genefull_features, 
            'GENE_SUMM': ss_gene_summary, 'GENEFULL_SUMM': ss_genefull_summary, 
            'BARCODE_STATS': ss_bc_stats, 'DEMUX': dem_file
            }
        
        for k, v in per_samp_check.items():
            if k == 'REG' and v == "":
                ss_dep_files = [
                    'REG', 'GENE_FEATURE', 'GENE_SUMM',
                    'GENEFULL_FEATURE', 'GENEFULL_SUMM',
                    'BARCODE_STATS'
                ]
                ss_dep_files = [ j for j in ss_dep_files if j not in samp_excl_progs ]
                samp_excl_progs.extend(ss_dep_files)

            elif k not in samp_excl_progs and v == "":
                samp_excl_progs.append(k)


        if len(samp_excl_progs) == len(per_samp_check):
            skip_sample=True

        if skip_sample:
            print(
                f"Skipping {sample} as all the input files "
                "are not present!!"
            )
            continue
           
        if isinstance(combo_log, list) or \
            (isinstance(combo_log, pd.DataFrame) and \
                not(combo_log['LAB']['SAMPLE']['SAMPLE'].str.contains(sample).any())
            ):

            # Add a kwargs style input for extra annotations
            row_list.append(write_logs(
                map_names.iloc[:, [2, 3, 1]].to_numpy().tolist(), map_names, files_dict, 
                samp_excl_progs, processed_data
                # sample=sample, 
                # prep=preparer, rep=replicate, 
                # set_num=set_val
                )
            )
        
        else:
            old_val = combo_log.loc[combo_log['LAB']['SAMPLE']['SAMPLE'].str.contains(sample), :].to_numpy().flatten().tolist()
            new_val = write_logs(
                        map_names.iloc[:, [2, 3, 1]].to_numpy().tolist(), map_names, files_dict, 
                        samp_excl_progs, processed_data
                        # sample=sample, 
                        # prep=preparer, rep=replicate, 
                        # set_num=set_val
                    )
            
            # If the new_val is different than the old value, replace it 
            if old_val != new_val:
                combo_log.loc[combo_log['LAB']['SAMPLE']['SAMPLE'].str.contains(sample)] = [new_val]
            else:
                print(
                    f"Skipping {sample} as the output file already "
                    "contains the same value for all columns!!"
                )
                continue

        print(f"Finished adding {sample} to the file")

    
    temp_df = pd.DataFrame(row_list, columns=pd.MultiIndex.from_frame(cols_list, names=["prog", "sub_prog", "curr_val"]))
    combo_log = pd.concat([combo_log.astype(temp_df.dtypes), temp_df.astype(combo_log.dtypes)], ignore_index=True)

    # Conversions for easier divisons
    combo_log[("STARsolo", "DEMUX", "N_CELLS_START")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX", "N_CELLS_START")]
        ]
    combo_log[("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT")]
        ]
    combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")]
        ]
    combo_log[("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS")]
        ]
    combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")]
        ]
    combo_log[("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS")]
        ]
    combo_log[("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS")]
        ]
    combo_log[("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS")]
        ]
    combo_log[("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS")]
        ]
    combo_log[("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS")] = [ 
        np.nan if i == '' or i == 'None' else float(i) for i in combo_log[("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS")]
        ]

    # If the demux file is not used for compilation then skip these steps
    if args.dem_info != None:
        try:
            # combo_log[("STARsolo", "DEMUX", "DOUBLET_PCT")] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_DOUBLET_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
            # combo_log[("STARsolo", "DEMUX_CS", "DOUBLET_PCT")] = combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")].astype(int)/combo_log[("STARsolo", "DEMUX", "N_CELLS_START")].astype(int)
            combo_log[("STARsolo", "DEMUX_CS", "DOUBLET_PCT")] = combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")].divide(combo_log[("STARsolo", "DEMUX", "N_CELLS_START")])
        except:
            print("Can't calculate Doublet ratio! Check output file for more info!")

        try:
            # combo_log[("STARsolo", "DEMUX", "NEGATIVE_PCT")] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_NEGATIVE_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
            combo_log[("STARsolo", "DEMUX_CS", "NEGATIVE_PCT")] = combo_log[("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS")].divide(combo_log[("STARsolo", "DEMUX", "N_CELLS_START")])
        except:
            print("Can't calculate Negative ratio! Check output file for more info!")

        try:
            combo_log[("STARsolo", "DEMUX_CS", "N_DEMUXED_CELLS")] = (combo_log[("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT")].subtract(combo_log[("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS")])
                        ).subtract(combo_log[("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS")])

        except:
            print("Can't calculate final demultiplexed cells from calico_solo run!")

        try:
            combo_log[("STARsolo", "DEMUX_CS", "CELL_RENTENTION")] = combo_log[("STARsolo", "DEMUX_CS", "N_DEMUXED_CELLS")].divide(combo_log[("STARsolo", "DEMUX", "N_CELLS_START")])
        except:
            print("Can't calculate percentage of cells retained after calico_solo demultiplexing")

        try:
            # combo_log[("STARsolo", "DEMUX", "DOUBLET_PCT")] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_DOUBLET_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
            combo_log[("STARsolo", "DEMUX_VS", "DOUBLET_PCT")] = combo_log[("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS")].divide(combo_log[("STARsolo", "DEMUX", "N_CELLS_START")])
        except:
            print("Can't calculate Doublet ratio! Check output file for more info!")

        try:
            # combo_log[("STARsolo", "DEMUX", "NEGATIVE_PCT")] = calc_ratio(combo_log["STARsolo"]["DEMUX"]["N_NEGATIVE_CELLS_CS"], combo_log["STARsolo"]["DEMUX"]["N_CELLS_START"])
            combo_log[("STARsolo", "DEMUX_VS", "NEGATIVE_PCT")] = combo_log[("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS")].divide(combo_log[("STARsolo", "DEMUX", "N_CELLS_START")])
        except:
            print("Can't calculate Negative ratio! Check output file for more info!")

        try:
            combo_log[("STARsolo", "DEMUX_VS", "N_DEMUXED_CELLS")] = (combo_log[("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT")].subtract(combo_log[("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS")])
                        ).subtract(combo_log[("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS")])

        except:
            print("Can't calculate final demultiplexed cells from vireo run")

        try:
            combo_log[("STARsolo", "DEMUX_VS", "CELL_RENTENTION")] = combo_log[("STARsolo", "DEMUX_VS", "N_DEMUXED_CELLS")].divide(combo_log[("STARsolo", "DEMUX", "N_CELLS_START")])
        except:
            print("Can't calculate percentage of cells retained after vireo demultiplexing")


    combo_log.replace([np.inf, -np.inf], np.nan, inplace=True)
    combo_log.to_csv(out, sep = "\t", index=False)


    sleep(30)
    

# Run this only when executed through Snakemake
if __name__ == "__main__":
    main()
