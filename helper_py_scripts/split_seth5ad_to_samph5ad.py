#!/usr/bin/env python
import re, glob2, os, argparse, gc
import pandas as pd, anndata as ad, scanpy as sc
from pathlib import Path
from collections import Counter

#parser = argparse.ArgumentParser(description="Split a(set of) h5ad(s) into multiple h5ad based on a obs column in it(them)")

#group1 = parser.add_argument_group("")
#parser.add_argument('inp_file', nargs='+', help="Paths(s) of h5ad files to split")
#parser.add_argument("-o", "--out_dir", default=Path("/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/h5ad"),help="Output directory (Directory will be created relative to current dir. Default path: /sc/arion/projects/psychAD/NPS-AD/freeze1_rc/h5ad)")
#parser.add_argument("--obs", default="SubID_cs", help="The obs column based on which h5ad file(s) will be split. Default obs name: SubID_cs")

#args = parser.parse_args()
'''
This script can't find changes in the constituents of the donor files but it can
identify an absence of a donor (either fresh or a newer application of demultiplexing methog
has resulted in the identification of a new donor)

'''

df = pd.read_csv('/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/New_wet_lab_spreadsheet.tsv', sep='\t')
out_dir = "/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/h5ad/"
for setid, sub in zip(df['set_ID'], df['SubID']):
	# If one of the donors is not present remove all donors and redo them
        # even if a set doesn't have a single good qualtiy donor this script will redo the whole set
	if not(os.path.isfile("{}{}_{}.h5ad".format(out_dir, setid[:-1], sub))) and glob2.glob("/sc/arion/projects/psychAD/final_count_matrix/solo/*_Sample-{}*_STARsolo_out.h5ad".format(setid[:-1])):
		samples = sorted(glob2.glob("/sc/arion/projects/psychAD/final_count_matrix/solo/*_Sample-{}*_STARsolo_out.h5ad".format(setid[:-1])))
		adatas = [ad.read(sample) for sample in samples] # create anndata array

		adata = ad.concat(adatas[:], index_unique='-', join='inner') # it will add a 'batch' variable with keys as its value and will inner join on the 'genes'
		output_files = adata.obs['SubID_cs'].unique()

		file_list = ["/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/h5ad/{}_{}.h5ad".format(setid[:-1], donor) for donor in output_files if donor not in ['Doublet', 'Negative', 'Not Present']]
		if all(map(os.path.isfile, file_list)):
			continue
		else:
                	map(os.remove, file_list)


		for op_f in output_files:
			if op_f not in ['Doublet', 'Negative', 'Not Present']:
				adata[adata.obs['SubID_cs'] == op_f].write("{}{}_{}.h5ad".format(out_dir, setid[:-1], op_f))

		del adata, adatas
	gc.collect()

	print("Processed Set: {}".format(setid[:-1]))
