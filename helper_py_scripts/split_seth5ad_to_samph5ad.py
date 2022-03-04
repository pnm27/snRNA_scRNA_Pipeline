#!/usr/bin/env python
import re, glob2, os, argparse, gc
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




df = pd.read_csv('/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/New_wet_lab_spreadsheet.tsv', sep='\t')
count_matrix_dir = "/sc/arion/projects/psychAD/final_count_matrix/solo/"
out_dir = "/sc/arion/projects/psychAD/NPS-AD/02_h5ad-by-donor/"

# Remove duplicated lines from the wet lab spreadsheet(remove HTO sets and replicate sets i.e. retain sets like A1-cDNA, C1-cDNA, S1-cDNA, etc)
df = df.loc[:, ['unique_sample_ID', 'SubID']]
df = df[df['unique_sample_ID'].str.endswith('cDNA')]
df.reset_index(inplace=True, drop=True)
df['set_ID'] = df['unique_sample_ID'].map(lambda x: x[:-6].replace('_', '-'))
df.drop_duplicates(subset=['set_ID'], ignore_index=True, inplace=True)

for setid, sub in zip(df['set_ID'], df['SubID']):

	all_files = glob2.glob(os.path.join(count_matrix_dir, "*_Sample-{}*_STARsolo_out.h5ad".format(setid)))
	all_donors = all([ os.path.isfile("{}{}_{}.h5ad".format(out_dir, setid, s)) for s in sub.split(',')])
	# If one of the donors is not present remove all donors and redo them
    # even if a set doesn't have a single good qualtiy donor this script will redo the whole set
    # Empty list is False
	if not all_donors and all_files:
		samples = sorted(all_files)
		adatas = [ad.read(sample) for sample in samples] # create anndata array

		adata = ad.concat(adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'
		output_files = adata.obs['SubID_cs'].unique()

		file_list = [os.path.join(out_dir, "{}_{}.h5ad".format(setid[:-6], donor)) for donor in output_files if donor not in ['Doublet', 'Negative', 'Not Present']]

		if all(map(os.path.isfile, file_list)) and check_same(adata, file_list):
			continue
		else:
			map(os.remove, file_list)


		for op_f in output_files:
			if op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" not in op_f:
				adata[adata.obs['SubID_cs'] == op_f].write("{}{}_{}.h5ad".format(out_dir, setid, op_f))
			elif op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" in op_f:
				new_name=op_f.replace('_unknown_donor', '_ud_')
				adata[adata.obs['SubID_cs'] == op_f].write("{}{}_{}.h5ad".format(out_dir, setid, new_name))

		del adata, adatas
		print("Processed incomplete Set: {}".format(setid))

	elif all_files and all_donors:
		samples = sorted(all_files)
		adatas = [ad.read(sample) for sample in samples] # create anndata array

		adata = ad.concat(adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'
		output_files = adata.obs['SubID_cs'].unique()

		file_list = [ os.path.join(out_dir, "{}_{}.h5ad".format(setid, s)) for s in sub.split(',')]

		old_adatas = [ad.read(f) for f in file_list] # create anndata array for the older per-donor h5ads
		old_adata = ad.concat(old_adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'

		del adatas, old_adatas
		gc.collect()
		if check_same(adata, ann_d2 = old_adata):
			print("Previously, completely processed set: {}".format(setid))
			del adata, old_adata
			gc.collect()
			continue
		else:
			map(os.remove, file_list)

		for op_f in output_files:
			if op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" not in op_f:
				adata[adata.obs['SubID_cs'] == op_f].write("{}{}_{}.h5ad".format(out_dir, setid, op_f))
			elif op_f not in ['Doublet', 'Negative', 'Not Present'] and "_unknown_donor" in op_f:
				new_name=op_f.replace('_unknown_donor', '_ud_')
				adata[adata.obs['SubID_cs'] == op_f].write("{}{}_{}.h5ad".format(out_dir, setid, new_name))

		del adata, old_adata
		print("Processed incomplete Set: {}".format(setid))

	else:
		print("Something looks wrong. PLease check either set number {} or donors {}".format(setid, sub))
		
	gc.collect()

	
