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
    if not all(sample + '-cDNA' in lf_df["LAB"]["SAMPLE"]["SAMPLE"] for sample in sample_name ):
        # print("The file {} doesn't contain all samples present in the \"fastq_files.txt\"!\nRemoving the file and producing a newer version of the file".format(log_file))
        os.remove(log_file)
        return log_file

    return []




# For splitting bams