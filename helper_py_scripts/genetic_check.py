#!/usr/bin/env python

import pandas as pd
from scipy.io import mmread
import vireoSNP
import numpy as np
import os
from collections import Counter
import warnings

def read_data(setID, round, demux_method):
	root = cellSNP_dir + '/Sample_' + setID + '-cDNA/'
	AD = mmread(root + 'cellSNP.tag.AD.mtx')
	DP = mmread(root + 'cellSNP.tag.DP.mtx')
	df_hto = pd.read_csv(demux_dir + demux_method + '/' + round + '_Sample-' + setID + '-cDNA.csv', index_col = 0)
	df_hto_cn = pd.read_csv(root + 'cellSNP.samples.tsv', sep = '\t', header = None)
	df_var = pd.read_csv(root + 'cellSNP.base.vcf.gz', sep = '\t', skiprows = 1)
	return AD, DP, df_hto, df_var

def slice_data(AD, DP, df_hto, hto, min_count = 20):
	ind_cell = df_hto.x.values == hto
	AD0 = AD.tocsc()[:, ind_cell]
	DP0 = DP.tocsc()[:, ind_cell]
	ind_fea = np.array(np.sum(DP0, axis = 1) >= min_count)[:, 0]
	AD0 = AD0.tocsr()[ind_fea]
	DP0 = DP0.tocsr()[ind_fea]
	df_var_hto = df_var.loc[ind_fea]
	return AD0, DP0, df_var_hto

def minmax(val_list):
	min_val = min(val_list)
	max_val = max(val_list)
	return (min_val, max_val)

def get_hto_var(ad_mat, dp_mat, df_hto, cell_var):
	hto_vars = {}
	cells = cell_var.columns[3:len(cell_var.columns)]
	#cells_v2 = []
	#for cell in cells:
		#cell2 = cell.replace('-1','')
		#cells_v2.append(cell)
	df_hto_v2 = df_hto[df_hto.index.isin(cells)]
	df_hto = df_hto_v2
	hto_list = np.unique(df_hto['x']).tolist()
	hto_list = [x for x in hto_list if "HTO" in x]
	for hto in hto_list:
		ad_mat0, dp_mat0, df_var_hto = slice_data(ad_mat, dp_mat, df_hto, hto, min_count = 10)
		id2var = [c for i, c in enumerate(zip(df_var_hto['#CHROM'], df_var_hto['POS']))]
		res = vireoSNP.vireo_wrap(ad_mat0, dp_mat0, n_donor=1, learn_GT=True,
								n_extra_donor=0, ASE_mode=False, fix_beta_sum=False,
								n_init=20, check_doublet=False, random_seed=None)
		id2var_df = pd.DataFrame(id2var, columns =['chr', 'position'])
		id2var_df['snp'] = id2var_df['chr'].astype(str) + '_' + id2var_df['position'].astype(str)
		id2var_df['snp'] = id2var_df['snp'].replace('chr','',regex = True)
		GT = []
		GT_prob = res['GT_prob']
		GT = np.repeat(np.nan, GT_prob.shape[0]).tolist()
		GT_prob = res['GT_prob']
		for i in range(0, GT_prob.shape[0]):
			for j in range(0, GT_prob.shape[2]):
				if GT_prob[i,0,j] > 0.5:
					GT[i] = j
					break
		id2var_df['GT'] = GT
		#remove heterozygote
		#id2var_df = id2var_df[id2var_df['GT'] != 1]
		id2var_df = id2var_df.dropna()
		id2var_df['GT'] = [int(x) for x in id2var_df['GT']]
		hto_vars[hto] = id2var_df
	return(hto_vars)

def same_var_pct(list1, list2):
	diflist = list1 - list2
	diflist = diflist.tolist()
	pct = diflist.count(0)/len(diflist)
	return(pct)

def best_match(hto, ref_var):
	df4check = ref_var[ref_var.snp.isin(hto_vars[hto]['snp'])]
	merged_df4check = pd.merge(left=hto_vars[hto], right=df4check, left_on='snp', right_on='snp')
	match_pct = []
	for i in range(4, len(merged_df4check.columns)):
		pct = same_var_pct(merged_df4check[merged_df4check.columns[3]], merged_df4check[merged_df4check.columns[i]])
		match_pct.append(pct)
	sampleID = merged_df4check.columns.tolist()[4:len(merged_df4check.columns)]
	match_pct_df = pd.DataFrame({'best_match_sampleID' : sampleID, 'match_pct' : match_pct})
	best_match_pct = match_pct_df[match_pct_df.match_pct == match_pct_df.match_pct.max()]
	best_match_pct.insert(0, 'hto', hto)
	best_match_pct['var_ct'] = merged_df4check.shape[0]
	return(best_match_pct)

def comp_df(set_ID, hto, ref_var, ref_source):
	df4check = ref_var[ref_var.snp.isin(hto_vars[hto]['snp'])]
	if df4check.shape[0] != 0:
		merged_df4check = pd.merge(left=hto_vars[hto], right=df4check, left_on='snp', right_on='snp')
		match_pct = []
		for i in range(4, len(merged_df4check.columns)):
			pct = same_var_pct(merged_df4check[merged_df4check.columns[3]], merged_df4check[merged_df4check.columns[i]])
			match_pct.append(pct)
		match_pct_boolean = ['FALSE']*len(match_pct)
		match_pct_boolean[match_pct.index(max(match_pct))] = 'TRUE'
		sampleID = merged_df4check.columns.tolist()[4:len(merged_df4check.columns)]
		set_ID_hto = set_ID + '_' + hto
		comp_match_pct_df = pd.DataFrame({'set_ID_hto' : set_ID_hto, 'SNParray_or_WGS_sampleID' : sampleID, 'match_pct' : match_pct, 'best_match' : match_pct_boolean, 'var_ct' : merged_df4check.shape[0], 'ref_source' : ref_source})
	else:
		comp_match_pct_df = pd.DataFrame({'set_ID_hto' : set_ID + hto, 'SNParray_or_WGS_sampleID' : ['NA'], 'match_pct' : ['NA'], 'best_match' : ['NA'], 'var_ct' : ['NA'], 'ref_source' : ref_source})
	return(comp_match_pct_df)

#reference_vars
ref_var = pd.read_csv('/sc/arion/projects/CommonMind/shan/CMC_TOPMED_IMPUTE/CMC_MSSM_NPSAD_genotyped_all.vcf.gz', sep = '\t', skiprows = 46)
ref_var_v2 = pd.DataFrame({'snp' : ref_var['#CHROM'].astype(str) + '_' + ref_var['POS'].astype(str)})
ref_var_v2['snp'] = ref_var_v2['snp'].replace("chr", "", regex=True)
for x in range(9, ref_var.shape[1]):
    ref_var_v2[ref_var.columns[x]] = ref_var[ref_var.columns[x]].replace(':.*','',regex=True)
ref_var_v3 = pd.DataFrame({'snp' : ref_var['#CHROM'].astype(str) + '_' + ref_var['POS'].astype(str)})
ref_var_v3['snp'] = ref_var_v3['snp'].replace("chr", "", regex=True)
for x in range(1, ref_var_v2.shape[1]):
    ref_var_v3[ref_var_v2.columns[x]] = np.where(ref_var_v2[ref_var_v2.columns[x]] == '0|0', 0, np.where(ref_var_v2[ref_var_v2.columns[x]] == '0|1', 1, np.where(ref_var_v2[ref_var_v2.columns[x]] == '1|0', 1, 2)))

round_num = 'round1'
demux_me = 'solo'
ref_source = 'CMC_SNParray'
cellSNP_dir = '/sc/arion/projects/psychAD/Single_cell_data/STARsolo_cellSNP/' + round_num + '/'
#set_IDs = os.listdir(cellSNP_dir)
#set_IDs = [x.replace('Sample_', '') for x in set_IDs]
#set_IDs = [x.replace('-cDNA', '') for x in set_IDs]
demux_dir = '/sc/arion/projects/psychAD/Single_cell_data/STARsolo_demux/'
demux_files = os.listdir(demux_dir + demux_me)
demux_files = [f for f in demux_files if round_num in f]
set_IDs = [x.replace(round_num + '_Sample-', '') for x in demux_files]
set_IDs = [x.replace('-cDNA.csv', '') for x in set_IDs]
genetic_check_dir = demux_dir + 'genetic_check/' + round_num + '/'
if not(os.path.isdir(genetic_check_dir)):
	os.mkdir(genetic_check_dir)
#demux_mes = ['htodemux', 'MULTIseq', 'solo']
comp_match_pct_df = pd.DataFrame(columns = ['set_ID_hto','SNParray_or_WGS_sampleID','match_pct','best_match','var_ct','ref_source'])
#set_ID = set_IDs[0]
for set_ID in set_IDs:
	ad_mat, dp_mat, df_hto, df_var = read_data(set_ID, round_num, demux_me)
	cell_var = pd.read_csv(cellSNP_dir + '/Sample_' + set_ID + '-cDNA/' + 'cellSNP.cells.vcf.gz', sep = '\t', skiprows = 37)
	cell_var = cell_var.drop(['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis = 1)
	cell_var.insert(2, 'snp', cell_var['#CHROM'].astype(str) + '_' + cell_var['POS'].astype(str))
	cell_var['snp'] = cell_var['snp'].replace('chr','',regex = True)
	#set_mapping_data = df4demux_comp[df4demux_comp['set_ID'] == df4demux_comp_uniq.at[i,'set_ID']]
	#set_mapping_data = set_mapping_data.reset_index()
	hto_vars = get_hto_var(ad_mat, dp_mat, df_hto, cell_var)
	#all comp df
	for hto in list(hto_vars.keys()):
		if hto_vars[hto].shape[0] != 0:
			comp_match_pct_df = comp_match_pct_df.append(comp_df(set_ID, hto, ref_var_v3, ref_source))
		else:
			df = pd.DataFrame({'set_ID_hto' : set_ID + hto, 'SNParray_or_WGS_sampleID' : ['NA'], 'match_pct' : ['NA'], 'best_match' : ['NA'], 'var_ct' : ['NA'], 'ref_source' : ref_source})
			comp_match_pct_df = comp_match_pct_df.append(df)
comp_match_pct_df.to_csv(genetic_check_dir + 'comp_match_pct_df.csv')

