#!/usr/bin/env python3

from __future__ import annotations
import pandas as pd
import os, json, logging, sys
import numpy as np
from time import sleep
import argparse, warnings # errno
from demultiplex_helper_funcs import (
    process_columns, auto_read, 
    has_wet_lab_value_column,
    get_demux_paths, process_swap_correction,
    get_filename, write_logs
)
from jsonschema import Draft202012Validator
from pathlib import Path


log_file = f"{Path(sys.argv[0]).stem}.log"
logger = logging.getLogger(__name__) # File only (DEFAULT)

# Prevent duplicate messages via the root logger
logger.propagate = True


def setup_logging(verbose: int):
    level = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO, 
        3: logging.DEBUG,
        }.get(verbose, logging.INFO)

    logging.basicConfig(
        level=level, filename=log_file,
        format="%(name)s - %(asctime)s - %(levelname)s - %(message)s",
        filemode='w',
    )


# This function is used to caluclate ratios like doublet pct and negative pct
def calc_ratio(numer, denom):
    ratios = []
    for i in range(len(numer)):
        ratios.append(int(numer[i])/int(denom[i]))
        
    return ratios


def get_all_columns(map_df: pd.DataFrame, 
    new_cols_to_add: list[list[str]]) -> pd.DataFrame:

    # Skip description (as present in the OG log files)
    # Reorder columns
    cols = map_df.iloc[:, [2, 3, 1]] 
    if new_cols_to_add:
        cl = pd.DataFrame(new_cols_to_add, columns=list(cols.columns.values))
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
        return cl
    else:
        return cols



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
    # Old deprecated optional argument
    parser.add_argument('-d', '--demul_dir', nargs='*',
        # nargs='?', 
        # help="Directory (comma-sep) containing demultiplexing stats. DEFAULT: current working directory",
        help=argparse.SUPPRESS, 
        # const=os.getcwd(), default=None
        )

    # Optional parameters
    parser.add_argument('-m', '--map_file', help="Mapping file that contains info on the headers in the output. DEFAULT: <current working directory>/Final_out_MAP_2.tsv", 
        default=os.path.join(os.getcwd(),"Final_out_MAP_2.tsv"))
    parser.add_argument('-o', '--output_file', help="output file. DEFAULT: <current working directory>/All_logs.tsv", default=os.path.join(os.getcwd(),"All_logs.tsv"))
    parser.add_argument('-b', '--bam_dir', help="Directory containing bam file(s). DEFAULT: current working directory", default=os.getcwd())
    parser.add_argument('-p', '--picard_dir', nargs='?', help="Directoy containing PICARD outputs. DEFAULT: current working directory", 
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
    parser.add_argument('--swap_correct', nargs='?', help="A csv or tsv file containing final swap corrected information. Will utilize info from here to figure which version of demultiplexing \
    was finally used for confirming donor assignments.", const="final_swap_correction.csv", default=None
    )
    parser.add_argument('-w', '--wet_lab_file', help="Path to file that contains either/all of: HTO info for each set, annotations, etc."
    )
    parser.add_argument('--verbose', '-v', action='count', default=0, 
        help="Increase output verbosity, default is ERROR and above "
        "(e.g. default for ERROR, -v for WARNING, -vv for INFO, -vvv for DEBUG)"
    )

    return parser



def main():
    """Main entry point"""

    # Parse arguments
    parser = get_argument_parser()
    args = parser.parse_args()
    
    setup_logging(args.verbose)

    # Parse all values
    bam_dir = args.bam_dir
    pic_dir = args.picard_dir
    dem_dir = args.demul_dir # DEPRACATED
    samples = args.samples if args.samples else args.samples_deprecated
    df = None # Default - Wet Lab File is not provided
    swap_corr_df = None # Default - Swap Correction File is not provided

    if dem_dir is not None:
        warnings.warn((
            "Positional argument for demultiplexing directory "
            "is deprecated and won't be used. Instead, specify in "
            "the file provided to key --common_annotations, now."), 
            DeprecationWarning
        )

    # Resolve which one to use
    if args.samples and args.samples_deprecated:
        parser.error("Provide either --samples or positional samples, not both.")

    if args.samples_deprecated:
        warnings.warn(
            "Positional 'samples' is deprecated. Use --samples instead.",
            DeprecationWarning
        )
  

    # If annotations is needed then provide the wet lab file
    if args.common_annotations is not None and args.wet_lab_file is None:
        warnings.warn(
            "Annotations are expected as json was provided but Wet lab file was NOT provided. \
            Annotations might be missing!",
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
        if k == 'demul_dir':
            continue
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

    except Exception as e:
        # Catch any other potential exceptions
        print(f"Will create a new file: {out}")

    if args.swap_correct is not None:
        swap_corr_df = auto_read(args.swap_correct)

    if args.common_annotations is not None:

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

        # Wet Lab file, Filter wet lab file's columns, if needed
        if has_wet_lab_value_column(data):
            df = auto_read(args.wet_lab_file)

    # Get all columns for the output df
    extra_cols = [ col["columnInLogs"] for col in data["columns"] if "columnInLogs" in col]
    cols_list = get_all_columns(map_names, extra_cols)
    demul_dirs = get_demux_paths(data)

    # Process each sample -----------------------------------------------------------------------------------------------------------------------------------------------
    # List containing per sample values as lists (list of lists)
    row_list = []
    for sample in samples:
        skip_sample = False # By default, don't skip any sample
        # NEW STYLE ANNOTIONS
        # PROCESS THE DATA
        processed_data = process_columns(data, sample, df).by_logs

        # create per sample copy of exclude_prog list
        samp_excl_progs = exclude_progs.copy()

        # Parse file structures for bam files, picard files and demultiplex info files
        bam_st = args.bam_struct.replace("<sample>", sample)
        pc_st = args.pc_struct.replace("<sample>", sample) if args.pc_struct is not None else ""
        dem_st = args.dem_struct.replace("<sample>", sample) if args.dem_struct is not None else ""
        

        dem_dir = process_swap_correction(data, swap_df=swap_corr_df,
                pool_name=sample, demux_paths=demul_dirs, logger=logger)
        
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
            "Demultiplex_stats": dem_file,
            # FUTURE
            # "Basic_Demultiplex_stats": dem_file_basic,
            # "Hashsolo_Demultiplex_stats": dem_file_hashsolo,
            # "Vireo_Demultiplex_stats": dem_file,
        }

        # Check for each sample in the list if it has all 
        # required files otherwise mark as "" for the respective sample (i.e. update samp_excl_progs list)
        per_samp_check = {
            'REG': ss_log_final, 'GC': pc_gc_file, 
            'RNASEQMETRIC': pc_rs_file, 'GENE_FEATURE': ss_gene_features, 
            'GENEFULL_FEATURE': ss_genefull_features, 
            'GENE_SUMM': ss_gene_summary, 'GENEFULL_SUMM': ss_genefull_summary, 
            'BARCODE_STATS': ss_bc_stats, 'DEMUX': dem_file
        }
        
        for k, v in per_samp_check.items():
            if k == 'REG' and not v:
                ss_dep_files = [
                    'REG', 'GENE_FEATURE', 'GENE_SUMM',
                    'GENEFULL_FEATURE', 'GENEFULL_SUMM',
                    'BARCODE_STATS'
                ]
                ss_dep_files = [ j for j in ss_dep_files if j not in samp_excl_progs ]
                samp_excl_progs.extend(ss_dep_files)

            elif k not in samp_excl_progs and not v:
                samp_excl_progs.append(k)


        if len(samp_excl_progs) == len(per_samp_check):
            skip_sample=True

        if skip_sample:
            print(
                f"Skipping {sample} as all the input files "
                "are not present!!"
            )
            continue

        row_list.append(write_logs(
            cols_list.to_numpy(), map_names, 
            files_dict, samp_excl_progs, sample, 
            logger=logger, 
            processed_data=processed_data,
            )
        )
    
    temp_df = pd.DataFrame(row_list, columns=pd.MultiIndex.from_frame(cols_list, names=["prog", "sub_prog", "curr_val"]))
    try:
        combo_log = pd.concat([combo_log.astype(temp_df.dtypes), temp_df.astype(combo_log.dtypes)], ignore_index=True, sort=False)
    except:
        combo_log = temp_df
    # resolve duplicates (keep latest)
    combo_log = combo_log.drop_duplicates(keep="last")

    # --- Convert columns to numeric ---
    cols_to_convert = [
        ("STARsolo", "DEMUX", "N_CELLS_START"),
        ("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT"),
        ("STARsolo", "DEMUX_CS", "N_DOUBLET_CELLS_CS"),
        ("STARsolo", "DEMUX_CS", "N_NEGATIVE_CELLS_CS"),
        ("STARsolo", "DEMUX_VS", "N_DOUBLET_CELLS_VS"),
        ("STARsolo", "DEMUX_VS", "N_NEGATIVE_CELLS_VS"),
    ]

    for col in cols_to_convert:
        combo_log[col] = pd.to_numeric(combo_log[col], errors="coerce")

    # If the demux file is not used for compilation then skip these steps
    if args.dem_info != None:
        base = combo_log[("STARsolo", "DEMUX", "N_CELLS_START")]
        mito = combo_log[("STARsolo", "DEMUX", "N_CELLS_LOW_MITO_PERCENT")]

        for mode in ["DEMUX_CS", "DEMUX_VS"]:
            try:
                doublets = combo_log[("STARsolo", mode, f"N_DOUBLET_CELLS_{mode.split('_')[-1]}")]
                negatives = combo_log[("STARsolo", mode, f"N_NEGATIVE_CELLS_{mode.split('_')[-1]}")]

                # Ratios
                combo_log[("STARsolo", mode, "DOUBLET_PCT")] = doublets.divide(base)
                combo_log[("STARsolo", mode, "NEGATIVE_PCT")] = negatives.divide(base)

                # Derived counts
                demuxed = mito.subtract(doublets).subtract(negatives)
                combo_log[("STARsolo", mode, "N_DEMUXED_CELLS")] = demuxed

                # Retention
                combo_log[("STARsolo", mode, "CELL_RENTENTION")] = demuxed.divide(base)

            except Exception as e:
                print(f"Error processing {mode}: {e}")


    combo_log.replace([np.inf, -np.inf], np.nan, inplace=True)
    combo_log.to_csv(out, sep = "\t", index=False)
    

# Run this only when executed through Snakemake
if __name__ == "__main__":
    main()
