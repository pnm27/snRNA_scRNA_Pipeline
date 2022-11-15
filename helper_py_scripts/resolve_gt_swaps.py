#!/sc/arion/work/prashf01/conda/envs/snakemake/bin/python

from typing import Union # Need verion > 3.5
import anndata as ad
import pandas as pd, os, gc

# assert sys.version_info >= (3, 5), "This script needs python version >= 3.5!"
# NPSAD-issue specific function ---------------------------------------------
def conv_type(x):
    try:
        return x.strip()
    except:
        return x

def get_prop_don(x, t_df) -> str:
    try:
        t_df.loc[t_df.iloc[:, 1:].isin([str(x).strip()]).any(axis=1), "SubID"].values[0]
        return x
    except:
        pass
# ---------------------------------------------------------------------------

# NPSAD-issue specific input-------------------------------------------------
conv_df = pd.read_csv("/sc/arion/projects/psychAD/pnm/NPSAD_gt_converter_Jaro_latest.csv", dtype=str)

all_cols=list(conv_df)
conv_df[all_cols] = conv_df[all_cols].applymap(conv_type)
# ---------------------------------------------------------------------------

ref_f=pd.read_csv("sample_to_donor_convertor_myUpdate_ver3.csv")
out_dir="/sc/arion/projects/psychAD/gt_swaps_resolved_donor_h5ad/"

for s in ref_f['SubID'].unique():
    # Only if the donor files are absent then create it
    # Modify it later
    if not os.path.isfile(os.path.join(out_dir, f"{s}.h5ad")):
        files_list_df = ref_f.loc[ref_f['SubID'] == s, ["SubID", "orig_pool", "file_n", "orig_id"]]
        adatas = []
        for i, j in files_list_df.iterrows():
            temp = ad.read(j['file_n'])
            if temp[temp.obs['SubID_vS'] == j['orig_id']].shape[0] != 0:
                temp = temp[temp.obs['SubID_vS'] == j['orig_id']].copy()
            else:
                temp_d = get_prop_don(temp.obs['SubID_vS'].unique(), conv_df)
                temp = temp[temp.obs['SubID_vS'] == temp_d].copy()

            temp.obs['SubID_vS'] = j['SubID']
            temp.obs['poolID_ref'] = j['orig_pool']
            adatas.append(temp)
            
        adata = ad.concat(adatas[:], index_unique='-', join='outer') # it will add a 'batch' variable with keys as its value and will outer join on the 'genes'
        
        adata.write(os.path.join(out_dir, f"{s}.h5ad"))

        del adatas, adata, temp
        gc.collect()