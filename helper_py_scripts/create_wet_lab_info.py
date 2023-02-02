#!/usr/bin/env python

#  The files_tracker file lets the user know if a particular wet lab file's (files') content has been changed at the level of (either):
# "unique_sample_ID", "hashtag", "ab_barcode"
# NOTE: After each run the files_tracker file will be overwritten
import os, re, sys, glob2, time
from openpyxl import load_workbook
from collections import Counter
import pandas as pd, numpy as np, argparse
from itertools import chain, repeat



def get_subid(ser1, df1, df2):
	ret_val=[]
	c=1
	wrong_sets=[]
	for set_n in ser1.unique():
		subids=[]
		# Using unique because of repeated sets
		for sample in df1.loc[df1["Set_number"] == set_n, "unique_sample_ID"].squeeze().unique():
			# print(df2.loc[df2["Sample_ID"] == sample, "SubID"].values[0])
			try:
				val=df2.loc[df2["Sample_ID"] == sample.strip(), "SubID"].values[0]
			except:
				val=str(set_n).replace('.0', '') + "_unknown_" + str(c)
				# wrong_sets.append((sample, set_n))
				c+=1
			subids.append(val)
	    
		ret_val.append(','.join(subids))

	# Assuming per set there will be 4 rows representing 2 reps plus cDNA and HTO
	return list(chain.from_iterable(repeat(v, 4) for v in ret_val))
        

def get_argument_parser():
	"""Generate and return argument parser."""
	#Parse Command-Line arguments
	parser = argparse.ArgumentParser(description="Demultiplex sample based on hahsolo produced output")

	parser.add_argument('input', help="(List of) input file(s). "
	"Assumed to be in xlsx format.", nargs='+',
	)
	parser.add_argument('-o', '--output', help="Output file name. "
	"Default: output.csv (in the current working dir).", 
	default="output.csv",
	)
	parser.add_argument('-c', '--converter', nargs='?', help="File hash created "
	"by Jaro (links 'set' and ''). Expected with 'tsv' extension. "
	"Default: output.tsv (in the current working dir).", 
	const="converter.xlsx", default=None
	)
	parser.add_argument('-l', '--files_tracker', help="This file retains "
	"info about the processed xlsx file(s). Expected with 'txt' extension, "
	"will overwrite the previous file. Default: files_tracker.txt "
	"(in the current working dir).", default="files_tracker.txt",
	)
	# parser.add_argument('--columns', nargs=5, help="List of column names "
	# "RESPECTIVELY to sample_ID (as present in the spreadsheets), HTO numbers, "
	# "SubIDs, and set number.",
	# metavar=('sample_ID', 'HTO_name', 'HTO_barcode', 'Sub_ID', 'Set_number'), 
	# default=['unique_sample_ID', 'hto', 'hto_barcode', 'SubID', 'Set_number'],
	# )
	parser.add_argument('-p', '--project_name', help="Project name with which "
	"'unique sample ID' column starts with",
	)

	return parser


def main():
	"""Main function"""

	parser = get_argument_parser()
	args = parser.parse_args()
	col_names = ["project", "sample_ID", "Rx", "Set_number", "unique_sample_ID",
			"sample_type", "sample_notes", "date_processed", "Preparer", "input",
			"DNA_stain", "Hash_pool_contents", "hashtag", "ab_barcode", "FACS_by",
			"cell_type", "yield", "FACS_notes", "10x_by", "post_spin_cell_conc",
			"perc_viable_10x", "dilution", "normalized_conc", "amount_in_10x", "targeted_recovery",
			"n_cells", "capture_method", "Bead_lot_num", "10x_notes"]

	# work_cols=args.columns


	try:
		op_df = pd.read_csv(args.output, sep='\t')

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
			t_df['unique_sample_ID'] = t_df['unique_sample_ID'].apply(lambda x: str(x) if isinstance(x, int) else x)

			# Strip leading and trailing whitespaces
			df_obj2 = t_df.select_dtypes(['object'])
			df_obj2 = df_obj2.applymap(str) # forcefully convert everything to string
			t_df[df_obj2.columns] = df_obj2.apply(lambda x: x.str.strip())
			
			# Save donor info
			t_donor_info = t_df.loc[~t_df["unique_sample_ID"].str.contains(args.project_name, na=False)]
			donor_info = pd.concat([donor_info, t_donor_info])

			# Temporarily add an empty column for SubID
			t_df["SubID"]=""
			# Retain multiplexed sample info
			t_df = t_df.loc[t_df["unique_sample_ID"].str.contains(args.project_name, na=False)]
			# if not t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"]).all():
			# For those samples that weren't present in the ouput already
			try:
				t_df = t_df[~ t_df["unique_sample_ID"].isin(op_df["unique_sample_ID"])]
			except:
				print(inp_f)
				print(op_df.columns)
				
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
				with open(args.files_tracker, 'a') as fout:
					fout.write("The file {} has been updated. The following samples inside this file were updated:".format(inp_f))
					for sample in df2["unique_sample_ID"]:
						fout.write("\t\t{}".format(sample))

		with open(args.files_tracker, 'a') as fout:
			fout.write("Hence, Re-do all related analyses.")
				

	else:
		t_df = pd.read_excel(args.input, names=col_names, skiprows=1, usecols=list(range(len(col_names))))
		t_df["filename"] = os.path.basename(args.input).replace('.xlsx', '')
		t_df['unique_sample_ID'] = t_df['unique_sample_ID'].apply(lambda x: str(x) if isinstance(x, int) else x)

		# Strip leading and trailing whitespaces
		df_obj2 = t_df.select_dtypes(['object'])
		df_obj2 = df_obj2.applymap(str) # forcefully convert everything to string
		t_df[df_obj2.columns] = df_obj2.apply(lambda x: x.str.strip())

		# Save donor info
		donor_info = t_df.loc[~t_df["unique_sample_ID"].str.contains(args.project_name, na=False)]
		# Temporarily add an empty column for SubID
		t_df["SubID"]=""
		
		# Retain multiplexed sample info
		t_df = t_df.loc[t_df["unique_sample_ID"].str.contains(args.project_name, na=False)]

		
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
			with open(args.files_tracker, 'a') as fout:
				fout.write("The file {} has been updated. The following samples inside this file were updated:".format(args.input))
				for sample in df2["unique_sample_ID"]:
					fout.write("\t\t{}".format(sample))
				fout.write("Hence, Re-do all related analyses.")


	if args.converter is not None:
		# Project dependent converter file
		conv_df = pd.read_excel(args.converter, usecols=list(range(15)))
		conv_df['Sample_ID'] = conv_df['Sample_ID'].apply(lambda x: str(x) if isinstance(x, int) else x)
		# Strip leading and trailing whitespaces
		df_obj = conv_df.select_dtypes(['object'])
		conv_df[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())


		op_df['SubID'] = get_subid(op_df['Set_number'], donor_info, conv_df)

	op_df['hashtag'] = op_df['hashtag'].apply(lambda x: re.sub('_', ',', x))
	op_df['ab_barcode'] = op_df['ab_barcode'].apply(lambda x: re.sub('_', ',', x))
	op_df['hashtag'] = op_df['hashtag'].apply(lambda x: x.replace('HTO#', ''))
	op_df.to_csv(args.output, sep='\t', index=False)

	ext = args.output[args.output.rfind('.'):]
	donor_filename = args.output[:args.output.rfind('.')] + '_donor' + ext
	donor_info['hashtag'] = donor_info['hashtag'].apply(lambda x: x.replace('#', ''))
	donor_info.to_csv(donor_filename, sep='\t', index=False)


if __name__ == '__main__':
    main()
