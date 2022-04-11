#!/usr/bin/env python
import re, glob2, os, argparse, gc, warnings
import pandas as pd, anndata as ad, scanpy as sc
from pathlib import Path
from collections import Counter


'''
This script can't find changes in the constituents of the donor files but it can
identify an absence of a donor (either fresh or a newer application of demultiplexing methog
has resulted in the identification of a new donor)

'''


# Redundant creation of h5ads. MAKE IT EFFECIENT!!
# Check whether 2 anndatas are same (per each donor too) or whether each set of per-donor file is same as the ones that will be produced by
def check_same(ann_d, files_l = None, ann_d2 = None) -> "boolean":
	if ann_d2 == None and files_l != None:
		d = [ os.path.basename(f).split('_')[1].replace('.h5ad', '').strip() for f in files_l ]
		for idx, vals in enumerate(d):
			t_ann = ad.read(files_l[idx])
			diff_by_shape = ann_d[ann_d.obs['SubID_cs'] == d].shape == t_ann.shape
			diff_by_reps =sorted(ann_d[ann_d.obs['SubID_cs'] == d].obs['rep'].unique()) == sorted(t_ann.obs['rep'].unique())
			return diff_by_reps and diff_by_reps


	elif ann_d2 != None and files_l == None:
		diff_by_shape = ann_d[(ann_d['SubID_cs'] != 'Doublet') & (ann_d.obs['SubID_cs'] != 'Negative') & (ann_d.obs['SubID_cs'] != 'Not Present')].shape == ann_d2.shape
		diff_by_reps =sorted(ann_d.obs['rep'].unique()) == sorted(ann_d2.obs['rep'].unique())
		diff_by_cell_comp = Counter(ann_d[(ann_d.obs['SubID_cs'] != 'Doublet') & (ann_d.obs['SubID_cs'] != 'Negative') & (ann_d.obs['SubID_cs'] != 'Not Present')].obs['SubID_cs']) == Counter(ann_d2.obs['SubID_cs'])
		return diff_by_shape and diff_by_reps and diff_by_cell_comp

	else:
		raise ValueError("Unrecognized combination for the function 'check_same'. Check parameters!")



# Read tsv, csv, txt (space-separated and tab-separated)
def read_files_ext(fname, lev) -> pd.DataFrame :
	if not os.path.isfile(fname):
		raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 

	if fname.endswith('.csv'):
		return pd.read_csv(fname, header=list(range(lev)))
	elif fname.endswith('.tsv'):
		return pd.read_csv(fname, sep='\t', header=list(range(lev)))
	elif fname.endswith('.txt'):
		warnings.warn("Guessing either space-sep or tab-sep can result in an undesired layout of the file")
		return pd.read_csv(fname, sep=r'\s', engine='python', header=list(range(lev)))
	else:
		raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")




def get_donors(inp_df, col_val1, col_val2, ds, mh, sn_f, sn_g, is_log) -> dict:
	dd={}
	if mh == 1: # Not a multiheader input file
		if ds == None: # Multiple lines exists for same pool value i.e. each donor of each pool is present in a unique row
			for pool in inp_df[col_val1].unique():
				temp_val = inp_df.loc[inp_df[col_val1] == pool, col_val2].tolist()
				if sn_f != None:
					key_name = sn_f.search(pool).group(sn_g) if sn_f.search(pool).group(sn_g) != None else pool
				else:
					key_name = pool

				dd[key_name] = list(map(str.strip, temp_val))

		elif len(ds) == 1: # One row (representing a pool/sample) contains all donor names separated by a SEP (ds)
			for pool in inp_df[col_val1]:
				temp_val = inp_df.loc[inp_df[col_val1] == pool, col_val2].values[0]
				if sn_f != None:
					key_name = sn_f.search(pool).group(sn_g) if sn_f.search(pool).group(sn_g) != None else pool
				else:
					key_name = pool

				# If the input file is the compiled log_file then donors are present as <donor_name>:<# of cell>
				if is_log:
					temp_val = [ str(v.split(':')[0]) for v in temp_val.split(ds) ]

				dd[key_name] = list(map(str.strip, temp_val))
		# elif ds == 'norm':
		# 	for sep in [' ', ',']
		# 		for pool in inp_df[col_val1]:
		# 			temp_val = inp_df.loc[inp_df[col_val1] == pool, col_val2].values[0]
		# 			dd[pool] = temp_val.split(ds)

		else:
			raise ValueError(f"Check the given \"donor_sep\" argument! Separator has length {len(ds)} ( expecting length 1)")

	else: # A multiheader input file
		if ds == None: # Multiple lines exists for same pool value i.e. each donor of each pool is present in a unique row
			for pool in inp_df[tuple(col_val1.split(","))].unique():
				temp_val = inp_df.loc[inp_df[tuple(col_val1.split(","))] == pool, tuple(col_val2.split(","))].tolist()
				if sn_f != None:
					key_name = sn_f.search(pool).group(sn_g) if sn_f.search(pool).group(sn_g) != None else pool
				else:
					key_name = pool

				dd[key_name] = list(map(str.strip, temp_val))

		elif len(ds) == 1: # One row (representing a pool/sample) contains all donor names separated by a SEP (ds)
			for pool in inp_df[tuple(col_val1.split(","))]:
				temp_val = inp_df.loc[inp_df[tuple(col_val1.split(","))] == pool, tuple(col_val2.split(","))].values[0]
				if sn_f != None:
					key_name = sn_f.search(pool).group(sn_g) if sn_f.search(pool).group(sn_g) != None else pool
				else:
					key_name = pool

				# If the input file is the compiled log_file then donors are present as <donor_name>:<# of cell>
				if is_log:
					temp_val = [ str(v.split(':')[0]) for v in temp_val.split(ds) ]

				dd[key_name] = list(map(str.strip, temp_val))

		else:
			raise ValueError(f"Check the given \"donor_sep\" argument! Separator has length {len(ds)} (expecting length 1)")

	return dd



# Parse arguments
parser = argparse.ArgumentParser(description="""\
Split each \"pooled\" h5ad into \"per-donor\" h5ad. Provide EITHER the compiled log file produced by update_logs.py script or the wet_lab_file provided to demul_samples.py script.
NOTE: If wet_lab_file is provided, certain related options (and same) provided to demul_samples.py script are required.
""", formatter_class=argparse.RawTextHelpFormatter)

# Required arguments
parser.add_argument('count_matrix_dir', help="""/
The directory containing all final count matrix files (h5ad). DEFAULT: \"<current_working_dir>/final_count_matrix/solo/\"
NOTE: Make sure the files are in the format (using the default value for this argument):
\"<current_working_dir>/final_count_matrix/solo/*<sample_name>*.h5ad\", where <sample_name> is the value in the column specified by the 3rd positional argument of this script
,if no optional argument \"sample_name_format\" is present otherwise <sample_name> is the value after formatting
""", default=os.path.join(os.getcwd(), "final_count_matrix/solo/"))
parser.add_argument('out_dir', help="The output directory for the donor-specific h5ad files. DEFAULT: \"<current_working_dir>/02_h5ad-by-donor/\"", default=os.path.join(os.getcwd(), "02_h5ad-by-donor/"))
parser.add_argument('sample_name_column', help="""\
Column header that contains sample/pool names in the input file (wet_lab_file or compiled log file). For multi/heirarchial headers specify the complete header as a single
string separating each level by comma.
""")


# Optional parameters
parser.add_argument('-l', '--log_file', nargs='?', help="Compiled log file produced by update_logs.py script.", const=os.path.join(os.getcwd()+"All_logs.tsv"), default=None)
parser.add_argument('-c', '--cols', nargs='?', help="""\
Column header that contains # of cells per donor for each pool. For multi/heirarchial headers (by default) specify the complete header as a single string separating each leve by comma.
DEFAULT:
  1) without any command-line value: \"STARsolo DEMUX N_CELLS_AFTER_DEMUX_CS\"
  2) without the argument: \"donor_names\"
  """, const="STARsolo, DEMUX, N_CELLS_AFTER_DEMUX_CS", default="donor_names")
parser.add_argument('--wet_lab_file', nargs='?', help="Path to file that contains HTO info for each set (either csv or tsv file)", const=os.path.join(os.getcwd()+ "metadata.csv"), default=None)
parser.add_argument('--multiheader', type=int, help="""
If the input file has heirarchial column headers specify last header line (assumed to be from the first line). EXAMPLE: If headers are present from first line to the third, 3.
DEFAULT: 1 (single level header)
""", default=1)
# parser.add_argument('--col_wl_file', nargs='?', help="Column header that contains # of cells per donor for each pool in the wet lab file. For multi/heirarchial headers specify the complete header as a single \
# 	string separating each leve by space. DEFAULT: \"donor\"", const="donor", default=None)
parser.add_argument('--donor_sep', help="If, per each sample in the input file, donors are all present in one row separated by some SEP then specify it here. Default: ' ' ", default=' ')
# parser.add_argument('--sample_name_format', help="""\
# Format/extract (regex with single group extracting the required name) for naming per donor h5ad files. Example: Files will be named as \"<out_dir>/<sample_name_format>_<donor_name>.h5ad\".
# DEFAULT: whole value present in the column specified.
# """, default=None)
# parser.add_argument('--samp_name_fmt_grp', type=int, help="If multiple groups are present in the \"sample_name_format\" argument then specify the number of the group. DEFAULT: No grouping", default=0)


args = parser.parse_args()

# Flag when both files are used as inputs
use_log=False

# sanity check
assert args.multiheader >= 1, "Value of the \"multiheader\" can not be lesser than 1 i.e. headers are expected in the input file"

# Check argument conditions
if args.log_file != None and args.wet_lab_file != None:
	warnings.warn("Compiled log file and wet lab file are both provided. Will only use compiled log file.")
	use_log=True
elif args.log_file == None and args.wet_lab_file == None:
	raise ValueError("Neither compiled log file nor wet_lab_file was provided. Please provide one of them!")
elif args.log_file != None and args.wet_lab_file == None and args.cols != None:
	use_log=True
	pass
elif args.log_file == None and args.wet_lab_file != None and args.cols != None:
	pass
else:
	raise ValueError("Given arguments are either insuficient or incomplete! Check the help!")


# Parse arguments
df = read_files_ext(args.log_file.strip(), args.multiheader) if use_log else read_files_ext(args.wet_lab_file.strip(), args.multiheader)
count_matrix_dir = args.count_matrix_dir
out_dir = args.out_dir
pool_regex = "-([0-9]+-[A-Za-z0-9]+)+"  # If don't want to extract then set to None
group_n = 1 # The group number to extract from regex in previous line

# Formatting pool/sample names
sn_fmt = re.compile(pool_regex) if pool_regex != None else None
sn_fmt_grps = group_n if sn_fmt != None else None


# If needed, process the input file
# Remove duplicated lines from the wet lab spreadsheet(remove HTO sets and replicate sets i.e. retain sets like A1-cDNA, C1-cDNA, S1-cDNA, etc)
# df = df.loc[:, ['unique_sample_ID', 'SubID']]
# df = df[df['unique_sample_ID'].str.endswith('cDNA')]
# df.reset_index(inplace=True, drop=True)
# df['set_ID'] = df['unique_sample_ID'].map(lambda x: x[:-6].replace('_', '-'))
# df.drop_duplicates(subset=['set_ID'], ignore_index=True, inplace=True)


# Check heirarchy levels for columns
if len(args.sample_name_column.strip().split(",")) == args.multiheader and len(args.cols.strip().split(",")) == args.multiheader: # Levels of columns should be same
	pass
elif len(args.sample_name_column.strip().split(",")) != args.multiheader and len(args.cols.strip().split(",")) == args.multiheader: # Expected Heirarchial column headers but invalid sample name columns' header
	raise ValueError(f"From the argument \'multiheader\' expected {args.multiheader} columns but third positional argument (sample_name column) is different: {len(args.sample_name_column.strip().split(','))}!")
elif len(args.sample_name_column.strip().split(",")) == args.multiheader and len(args.cols.strip().split(",")) != args.multiheader: # Expected Heirarchial column headers but invalid donor name columns' header
	raise ValueError(f"From the argument \'multiheader\' expected {args.multiheader} columns but the argument \'cols\' has different levels: {len(args.cols.strip().split(','))}!")
else: # Expected Heirarchial column headers but invalid donor name columns' header
	raise ValueError(f"The argument \'multiheader\' has a value {args.multiheader} but the third positional argument (sample_name column) and argument \'cols\' have different levels: \
		{len(args.sample_name_column.strip().split(','))} and {len(args.cols.strip().split(','))}, respectively!")



# Check whether all columns exist in the given df or not
# Check if pool/sample name columns isn't duplicated
if args.multiheader > 1:
	assert tuple(args.sample_name_column.strip().split(",")) in df.columns, "Sample_name column is absent in the provided input file! Check heirarchy order or individual headers"
	assert all(~df.duplicated(subset=tuple(args.sample_name_column.strip().split(",")))), "Sample_name_column contains a duplicated pool/sample value!"
	assert tuple(args.cols.strip().split(",")) in df.columns, "Donor column is absent in the provided input file! Check heirarchy order or individual headers"
else:
	assert args.sample_name_column.strip() in df.columns, "Sample_name column is absent in the provided input file!"
	assert all(~df.duplicated(subset=args.sample_name_column.strip().split())), "Sample_name_column contains a duplicated pool/sample value!"
	assert args.cols.strip() in df.columns, "Donor column is absent in the provided input file!"	


# Check if sample_formatting returns a value or not, if present
if sn_fmt != None and args.multiheader > 1:
	assert all(df[tuple(args.sample_name_column.strip().split(","))].apply(lambda x: sn_fmt.search(x).group(sn_fmt_grps)) != None), "Sample name formatting seems to not return any value! Check the arguments \
	\'sample_name_format\', \'sample_name_column\' and \'samp_name_fmt_grp\'"


# Get donors as dict (with sample_names as keys)
donors_dict = get_donors(df, args.sample_name_column, args.cols, args.donor_sep, args.multiheader, sn_fmt, sn_fmt_grps, use_log)


# Donor names are according to the annotation present in the "obs" column of the respective 'pooled' anndata (final count matrix)
for p, d in donors_dict.items():
	# If reps are present specify how each rep differs like:
	# all_files = glob2.glob(os.path.join(count_matrix_dir, "*_Sample-{}*_STARsolo_out.h5ad".format(setid)))\
	all_files = glob2.glob(os.path.join(count_matrix_dir, f"*{p}*/*.h5ad"))

	all_donors = all([ os.path.isfile(os.path.join(out_dir, f"{p}_{a}.h5ad")) for a in d]) if isinstance(d, list) else os.path.isfile(os.path.join(out_dir, f"{p}_{d}.h5ad"))

	# If one of the donors is not present remove all donors and redo them
    # even if a set doesn't have a single good qualtiy donor this script will redo the whole set
    # Empty list is False
	if not all_donors and all_files:
		samples = sorted(all_files)
		adatas = [ad.read(sample) for sample in samples] # create anndata array

		adata = ad.concat(adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'
		output_files = adata.obs['SubID_cs'].unique()

		file_list = [os.path.join(out_dir, f"{p}_{donor}.h5ad") for donor in output_files if donor not in ['Doublet', 'Negative', 'Not Present']]

		# For 2 situations:
		# 1) the donors present in the final count matrix is lesser than expected i.e. if pools lost donors after demultiplexing
		# 2) exactly "same" per-donor h5ad files are created
		if all(map(os.path.isfile, file_list)) and check_same(adata, file_list):
			continue
		else:
			map(os.remove, file_list)


		for op_f in output_files:
			if op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" not in op_f:
				adata[adata.obs['SubID_cs'] == op_f].write(os.path.join(out_dir, "{p}_{op_f}.h5ad"))
			elif op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" in op_f:
				new_name=op_f.replace('_unknown_donor', '_ud_')
				adata[adata.obs['SubID_cs'] == op_f].write(os.path.join(out_dir, "{p}_{op_f}.h5ad"))

		del adata, adatas
		print("Processed incomplete Set: {}".format(setid))

	elif all_files and all_donors: # Check when everything is unchanged i.e. all input files are present and all per-donor h5ad files are present too
		samples = sorted(all_files)
		adatas = [ad.read(sample) for sample in samples] # create anndata array

		adata = ad.concat(adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'
		output_files = adata.obs['SubID_cs'].unique()

		file_list = [ os.path.join(out_dir, f"{p}_{donor}.h5ad") for donor in output_files if donor not in ['Doublet', 'Negative', 'Not Present']]

		old_adatas = [ad.read(f) for f in file_list] # create anndata array for the older per-donor h5ads
		old_adata = ad.concat(old_adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'

		del adatas, old_adatas
		gc.collect()
		if check_same(adata, ann_d2 = old_adata):
			print(f"Previously, completely processed set: {p}")
			del adata, old_adata
			gc.collect()
			continue
		else:
			map(os.remove, file_list)

		for op_f in output_files:
			if op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" not in op_f:
				adata[adata.obs['SubID_cs'] == op_f].write(os.path.join(out_dir, "{p}_{op_f}.h5ad"))
			elif op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" in op_f:
				new_name=op_f.replace('_unknown_donor', '_ud_')
				adata[adata.obs['SubID_cs'] == op_f].write(os.path.join(out_dir, "{p}_{op_f}.h5ad"))

		del adata, old_adata
		print(f"Processed incomplete Set: {p}")

	else:
		print("Something looks wrong. PLease check either set number {} or donors {}".format(p, d))
		
	gc.collect()

