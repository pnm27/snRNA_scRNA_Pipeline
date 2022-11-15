#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

from typing import Union # Need verion > 3.5
import anndata as ad
import pandas as pd, os

# assert sys.version_info >= (3, 5), "This script needs python version >= 3.5!"


ref_f=pd.read_csv("sample_to_donor_convertor_myUpdate_ver3.csv")
out_dir="/sc/arion/projects/psychAD/gt_swaps_resolved_donor_h5ad/"

for s in ref_f['SubID'].unique():
    files_list_df = ref_f.loc[ref_f['SubID'] == s, ["SubID", "orig_pool", "file_n", "orig_id"]]
    adatas = []
    for i, j in files_list_df.iterrows():
        temp = ad.read(j['file_n'])
        temp = temp[temp.obs['SubID_vS'] == j['orig_id']].copy()
        temp.obs['SubID_vS'] = j['SubID']
        temp.obs['poolID_ref'] = j['orig_pool']
        adatas.append(temp)
        
    adata = ad.concat(adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'
    
    adata.write(os.path.join(out_dir, f"{s}.h5ad"))

    del adatas, adata, temp
    gc.collect()