#!/usr/bin/env python

#  The files_tracker file lets the user know if a particular wet lab file's (files') content has been changed at the level of (either):
# "unique_sample_ID", "hashtag", "ab_barcode"
# NOTE: After each run the files_tracker file will be overwritten
import os, re, sys, glob2, time
from openpyxl import load_workbook
from collections import Counter
import pandas as pd, numpy as np, argparse



def get_subid(ser1, df1, df2):
	ret_val=[]
	c=1
	wrong_sets=[]
	for set_n in ser1:
	    subids=[]
        # Using unique because of repeated sets
	    for sample in df1.loc[df1["Set_number"] == set_n, "unique_sample_ID"].squeeze().unique():
	#             print(df2.loc[df2["Sample_ID"] == sample, "SubID"].values[0])
	        try:
	            val=df2.loc[df2["Sample_ID"] == sample, "SubID"].values[0]
	        except:
	            val=str(set_n).replace('.0', '') + "_unknown_donor" + str(c)
	            # wrong_sets.append((sample, set_n))
	            c+=1
	        subids.append(val)
	    
	    ret_val.append(','.join(subids))

	#     print(c)
	# print(set(wrong_sets))
	return ret_val
        


#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Demultiplex sample based on hahsolo produced output")

parser.add_argument('input', help="(List of) input file(s). Assumed to be in xlsx format.", nargs='+')
parser.add_argument('-o', '--output', help="Output file name. Default: output.csv (in the current working dir).", default="output.csv")
parser.add_argument('-c', '--converter', help="File hash created by Jaro (links 'set' and ''). Expected with 'tsv' extension. Default: output.tsv (in the current working dir).", default="output.csv")
parser.add_argument('-l', '--files_tracker', help="This file retains info about the processed xlsx file(s). Expected with 'txt' extension, will overwrite the previous file. Default: files_tracker.txt (in the current working dir).", default="files_tracker.txt")
parser.add_argument('--columns', nargs=5, help="List of column names RESPECTIVELY to sample_ID (as present in the spreadsheets), \
    HTO numbers, SubIDs, and set number.", metavar=('sample_ID', 'HTO_name', 'HTO_barcode', 'Sub_ID', 'Set_number'), default=['unique_sample_ID', 'hto', 'hto_barcode', 'SubID', 'Set_number'])

args = parser.parse_args()
col_names = ["project", "sample_ID", "Rx", "Set_number", "unique_sample_ID",
	       "sample_type", "sample_notes", "date_processed", "Preparer", "input",
           "DNA_stain", "Hash_pool_contents", "hashtag", "ab_barcode", "FACS_by",
           "cell_type", "yield", "FACS_notes", "10x_by", "post_spin_cell_conc",
           "perc_viable_10x", "dilution", "normalized_conc", "amount_in_10x", "targeted_recovery",
           "n_cells", "capture_method", "Bead_lot_num", "10x_notes"]

work_cols=args.columns


try:
	op_df = pd.read_csv(args.output)

except:
	op_df = pd.DataFrame(columns=col_names)
	op_df["filename"]=""
	op_df["SubID"]=""

# Handling sinlge-file as an input or multiple files as input
if isinstance(args.input, list):
	donor_info=pd.DataFrame(columns=col_names)
	donor_info["filename"]=""
	for inp_f in args.input:
		t_df = pd.read_excel(inp_f, names=col_names, skiprows=1, usecols=list(range(len(col_names))))
		t_df["filename"] = os.path.basename(inp_f).replace('.xlsx', '')
		# Save donor info
		t_donor_info = t_df.loc[~t_df["unique_sample_ID"].str.contains("NPSAD", na=False)]
		donor_info = pd.concat([donor_info, t_donor_info])
		# Temporarily add an empty column for SubID
		t_df["SubID"]=""
		# Retain multiplexed sample info
		t_df = t_df.loc[t_df["unique_sample_ID"].str.contains("NPSAD", na=False)]
		# if not t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"]).all():
		# For those samples that weren't present in the ouput already
		t_df = t_df[~ t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"])]
		if not t_df.empty:
			op_df = pd.concat([op_df, t_df])
		# For those samples that were present in the output
		t_df_present = t_df[t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"])]
		# If the previous df is a subset of the output i.e. these samples were exactly the same previously
		# at the level specified at beginning of this script
		df1=op_df.loc[op_df["unique_sample_ID"].isin(t_df["unique_sample_ID"]), ["unique_sample_ID", "hashtag", "ab_barcode"]]
		df2=t_df_present[["unique_sample_ID", "hashtag", "ab_barcode"]]
		if df1.merge(df2).shape != df1.shape:
			op_df.drop(op_df[op_df["unique_sample_ID"].isin(t_df["unique_sample_ID"])].index, inplace=True)
			op_df = pd.concat([op_df, t_df_present])
			with open(args.files_tracker, 'w') as fout:
				fout.write("The file {} has been updated. The following samples inside this file were updated:".format(inp_f))
				for sample in df2["unique_sample_ID"]:
					fout.write("\t\t{}".format(sample))

	with open(args.files_tracker, 'w') as fout:
		fout.write("Hence, Re-do all related analyses.")
			

else:
	t_df = pd.read_excel(args.input, names=col_names, skiprows=1, usecols=list(range(len(col_names))))
	t_df["filename"] = os.path.basename(args.input).replace('.xlsx', '')
	# Save donor info
	donor_info = t_df.loc[~t_df["unique_sample_ID"].str.contains("NPSAD", na=False)]
	# Temporarily add an empty column for SubID
	t_df["SubID"]=""
	# Retain multiplexed sample info
	t_df = t_df.loc[t_df["unique_sample_ID"].str.contains("NPSAD", na=False)]
	# if not t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"]).all():
	# For those samples that weren't present in the ouput already
	t_df = t_df[~ t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"])]
	if not t_df.empty:
		op_df = pd.concat([op_df, t_df])
	# For those samples that were present in the output
	t_df_present = t_df[t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"])]
	# If the previous df is a subset of the output i.e. these samples were exactly the same previously
	# at the level specified at beginning of this script
	df1=op_df.loc[op_df["unique_sample_ID"].isin(t_df["unique_sample_ID"]), ["unique_sample_ID", "hashtag", "ab_barcode"]]
	df2=t_df_present[["unique_sample_ID", "hashtag", "ab_barcode"]]
	if df1.merge(df2).shape != df1.shape:
		op_df.drop(op_df[op_df["unique_sample_ID"].isin(t_df["unique_sample_ID"])].index, inplace=True)
		op_df = pd.concat([op_df, t_df_present])
		with open(args.files_tracker, 'w') as fout:
			fout.write("The file {} has been updated. The following samples inside this file were updated:".format(args.input))
			for sample in df2["unique_sample_ID"]:
				fout.write("\t\t{}".format(sample))
			fout.write("Hence, Re-do all related analyses.")



conv_df = pd.read_excel(args.converter, usecols=list(range(15)))
# conv_df.rename({'xxx':'Set_number'}, axis=1, inplace=True)

# n_df = op_df.loc[op_df["unique_sample_ID"].str.contains("NPSAD", na=False)]
# n_df = n_df.drop_duplicates(ignore_index=True)
# n_df = n_df.reset_index(drop=True)

op_df['SubID'] = get_subid(op_df['Set_number'], donor_info, conv_df)

op_df['hashtag'] = op_df['hashtag'].apply(lambda x: re.sub('_', ',', x))
op_df['ab_barcode'] = op_df['ab_barcode'].apply(lambda x: re.sub('_', ',', x))
# n_df = n_df.merge(conv_df, on='Set_number')

# n_df2 = n_df.drop(list(range(152, 156))+list(range(192, 196)))
# n_df2 = n_df2.reset_index(drop=True)


# n_df2.to_csv(args.output, sep='\t', index=False)
op_df.to_csv(args.output, sep='\t', index=False)