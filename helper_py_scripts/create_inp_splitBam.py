#!/usr/bin/env python
from itertools import repeat
import pandas as pd, numpy as np, anndata as ad, scanpy as sc
import sys, os, re, glob2, argparse
from collections import Counter


# Save dataframe to file
def save_df(df, suff, op):
	if suff.lower() == '.csv':
		df.to_csv(op, index=False)
	elif suff.lower() == '.tsv' or suff.lower() == '.txt':
		df.to_csv(op, sep = "\t", index=False)
	else:
		print("Can't save to file because of incorrect extension. Extension can be 'csv', 'tsv' or 'txt'")


def get_argument_parser():
	"""Generate and return argument parser."""

	#Parse Command-Line arguments
	parser = argparse.ArgumentParser(description="Create a two-columned txt file with barcodes corresponding to each Donor when provided with annotated h5ad")

	parser.add_argument('inp', help="Path to cached final count matrix with annotated barcodes")
	parser.add_argument('output', help="Provide output file name.")

	# Optional parameters

	parser.add_argument('--overwrite', action='store_true', help="Flag describing whether to overwrite or not. If not, \"_2\" will be appended")

	return parser


def main():


	parser = get_argument_parser()
	# Parse arguments
	args = parser.parse_args()

	inp_h5ad = args.inp
	fout = args.output
	rem_op = args.overwrite


	file_ext = re.search(r'(\.[^.]+)$', fout).group(1)

	adata = ad.read(inp_h5ad)

	d = Counter(adata.obs['SubID_cs'])
	for vals in ['Doublet', 'Negative', 'Not Present']:
		d.pop(vals, None)

	df_l = []
	for k in d.keys():
		df_l.extend( (k, bc) for bc in adata[adata.obs['SubID_cs'] == k].obs_names.to_list() )

	temp_df = pd.DataFrame(df_l, columns=['Subj_ID', 'barcodes'])

	if rem_op:
		try:
			os.remove(fout)

		except:
			pass
		save_df(temp_df, file_ext, fout)

	else:	
		file_pref = re.search(r'(.*)\.[^.]+$', fout).group(1)
		new_name = file_pref + '_2' + file_ext
		save_df(temp_df, file_ext, new_name)


if __name__ == '__main__':
	main()