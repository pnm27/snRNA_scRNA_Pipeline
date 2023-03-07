#!/usr/bin/env python3

""" Creates
"""

from itertools import repeat
import pandas as pd, argparse, os, re
from collections import Counter


# Save dataframe to file
def save_df(df, suff, op, sep=None):
	if suff.lower() == '.csv':
		df.to_csv(op, index=False)
	elif suff.lower() == '.tsv':
		df.to_csv(op, sep = "\t", index=False)
	elif suff.lower() == '.txt':
		if sep is not None:
			df.to_csv(op, sep = sep, index=False)
		else:
			df.to_csv(op, sep = " ", index=False)
	else:
		print(
			"Can't save to file because of incorrect extension. "
			"Extension can be 'csv', 'tsv' or 'txt'."
			)


def get_argument_parser():
	"""Generate and return argument parser."""
    
	# Parse Command-Line arguments
	parser = argparse.ArgumentParser(description="Produce inputs (text files) per sample for cellSNP"
	", containing filtered barcodes for each. If 'prev' flag is used then an h5ad file is expected"
	" and if not then a path containing the 10 mtx files. For an h5ad file one can specifiy if "
	"cell classifications are present or not and if present, then whether certain cell classes "
	"need to be removed or not ('keep_all_cells' flag)")

	parser.add_argument('vireo_inp', help="Path to donor_ids.tsv "
	"matrix (h5ad) or Path containing 10x mtx files.",
	)
	parser.add_argument('output', help="Path to store the barcodes "
    "file"
    )
	parser.add_argument('--converter_file', help="If names from vireo output "
    "needs to be changed"
    )
    
	return parser



def main():

	parser = get_argument_parser()
	args = parser.parse_args()

	# Store output_file and Create necessary folders
	op = args.output

	# Create parent dir(s) to the output
	if not os.path.isdir(op.replace('/' + os.path.basename(op), '')):
		os.makedirs(op.replace(os.path.basename(op), ''))
	
	vir_class = pd.read_csv(args.vireo_inp, sep='\t', usecols=["cell", "donor_id"])
	vir_class = vir_class[(vir_class["donor_id"] != "doublet") & (vir_class["donor_id"] != "unassigned")]
	vir_class.reset_index(drop=True, inplace=True)
	# vir_class.rename(columns={"cell":"barcodes", "donor_id":"Subj_ID"}, inplace=True, errors="raise")

	# For converting genotype IDs to Subject IDs-------------------------------
	if args.converter_file is not None:
		# Project-specific filtering of donor_ids for conversion
		vir_class["donor_id"] = vir_class["donor_id"].apply(lambda x: '_'.join(x.split('_')[1:]) if x.startswith('0_') else x)

		conv_df = pd.read_csv(args.converter_file)
		conv_df = conv_df.loc[conv_df['primary_genotype2'].isin(vir_class['donor_id']), ["SubID", "primary_genotype2"]]
		vir_class['Subj_ID'] = vir_class['donor_id'].apply(lambda x: conv_df.loc[conv_df["primary_genotype2"] == x, "SubID"].values[0] if not x.startswith('donor') else x)
		del vir_class['donor_id']
		vir_class.rename(columns={"cell":"barcodes"}, inplace=True, errors="raise")
		
	# -------------------------------------------------------------------------
	else:
		vir_class.rename(columns={"cell":"barcodes", "donor_id":"Subj_ID"}, inplace=True, errors="raise")

	vir_class = vir_class[["Subj_ID", "barcodes"]]

	file_ext = re.search(r'(\.[^.]+)$', op).group(1)

	save_df(vir_class, file_ext, op)


if __name__ == '__main__':

    main()
    sleep(60)
    