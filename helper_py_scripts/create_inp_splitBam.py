#!/usr/bin/env python3
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
		print("Can't save to file because of incorrect extension. "
	"Extension can be 'csv', 'tsv' or 'txt'")


def get_argument_parser():
	"""Generate and return argument parser."""

	#Parse Command-Line arguments
	parser = argparse.ArgumentParser(
		description="Create a two-columned txt file with barcodes "
		"corresponding to each Donor when provided with annotated h5ad."
		)
	
	parser.add_argument('inp', help="Path to cached final count matrix "
		     "with annotated barcodes")
	parser.add_argument('output', help="Provide output file name with "
		     "extension. Supports only csv, tsv or txt (tab sep).")

	# Optional parameters

	parser.add_argument('--overwrite', action='store_true', 
			help="Flag describing whether to overwrite or not. "
			"If not, \"_2\" will be appended."
			)
	
	parser.add_argument('--converter', 
			help="If names in the h5ad classification need to be changed."
			" Only csv file is supported."
			)

	parser.add_argument('--split_by', nargs='+', 
			help="To split bams by classifications preoduced by either "
			"one or both of calico_solo and vireo provide column(s) "
			"that contain the classification. If vireo, specify 'vireo'.",
			metavar="method", dest="method",
			)
	
	parser.add_argument('--demux_suffix', nargs='+', 
			help="Suffix(es) to use for naming the output files according "
			"to the demultiplexing methods specified.",
			metavar="suffix",
			)
	

	return parser


def main():
	r"""Main entry function

	This function creates a 2-columned text file with donors as the 
	first column with its corresponding barcodes in the second.
	If multiple demultiplexing softwares have annotated the cells
	and one wants to create sep output files for each of them
	"""

	parser = get_argument_parser()
	# Parse arguments
	args = parser.parse_args()

	assert len(args.method) == len(args.demux_suffix), "Number of suffixes"
	" and number of demux methods aren't same!"

	inp_h5ad = args.inp
	fout = args.output
	rem_op = args.overwrite
	suff = [ s[1:] if s.startswith('_') else s for s in args.demux_suffix]


	file_ext = re.search(r'(\.[^.]+)$', fout).group(1)
	filename_prefix = fout.replace(file_ext, '')

	adata = ad.read(inp_h5ad)

	for i, dem in enumerate(args.method):

		d = Counter(adata.obs[dem])

		if "vs" in dem or "vireo" in dem:
			for vals in ['Doublet', 'Negative', 'Not Present']:
				d.pop(vals, None)
		else:
			for vals in ['doublet', 'unassigned']:
				d.pop(vals, None)

		df_l = []
		for k in d.keys():
			df_l.extend( (k, bc) for bc in adata[adata.obs[dem] == k].obs_names.to_list() )

		temp_df = pd.DataFrame(df_l, columns=['Subj_ID', 'barcodes'])
		op_f = filename_prefix + '_' + suff[i] + file_ext

		if rem_op:
			try:
				os.remove(op_f)

			except:
				pass
			save_df(temp_df, file_ext, op_f)

		else:	
			file_pref = filename_prefix + '_' + suff[i] + '_2'
			new_name = file_pref + '_2' + file_ext
			save_df(temp_df, file_ext, new_name)


if __name__ == '__main__':
	main()