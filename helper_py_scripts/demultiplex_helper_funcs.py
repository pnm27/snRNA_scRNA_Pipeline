#!/usr/bin/env python3

import pandas as pd, numpy as np, os, sys
from collections import OrderedDict as ord_dict
from typing import Union, Optional, Literal # Need verion > 3.5


assert sys.version_info >= (3, 5), "This script needs python version >= 3.5!"

# Basic helper functions-----------------------------------------------------
def auto_read(fname, lev=1, **kwargs) -> pd.DataFrame :
    if not os.path.isfile(fname):
        raise OSError(
            f"The given file {fname} doesn't exist and annotations are "
            "impossible without this file!"
            ) 
    if fname.endswith('.csv'):
        if lev == 1:
            return pd.read_csv(fname, header=0, **kwargs)
        else:
            return pd.read_csv(fname, header=list(range(lev)), **kwargs)
    elif fname.endswith('.tsv') or fname.endswith('.txt'):
        if lev == 1:
            return pd.read_csv(fname, sep='\t', header=0, **kwargs)
        else:
            return pd.read_csv(fname, sep='\t', header=list(range(lev)), **kwargs)
    else:
        raise OSError(
            f"The given file {fname} doen't have either csv or tsv "
            "extension. Other extensions are not supported!"
            )


# Simplify these 'parse' functions
def parse_file(wet_lab_df, cols, s_name, hs, d_con) \
    -> Union[list[str], tuple[list[str], list[str]]]:

    hto_col = cols[0]
    subid_col = cols[1]

    ret_list = []
    if isinstance(hto_col, str):
        htos = wet_lab_df[hto_col]
        test_len = len(htos)
        if not d_con:
            sub = wet_lab_df[subid_col]
            assert len(htos) == len(sub), "SubID and HTO length not equal"
        # test_len2 = len(sub)

        # No command-line params and inference that all HTOs are present 
        # in one row separated by ","
        if hs == None and test_len == 1 and \
            htos.str.count(',').values[0] > 1:
            if not d_con:
                assert sub.str.count(',').values[0] == (
                    htos.str.count(',').values[0]
                ), "Not equal number of hto and donor names"
                ret_list.append(htos.values[0].split(','))
                ret_list.append(sub.values[0].split(','))
            else:
                ret_list.extend(htos.values[0].split(','))

            return ret_list
        # No command-line params and inference that all HTOs are present 
        # in one row separated by whitespaces    
        elif hs == None and test_len == 1 and \
            len(htos.split()) > 1:
            
            if not d_con:
                assert sub.str.count().values[0] == (
                    htos.str.count().values[0]
                ), "Not equal number of hto and donor names"
                ret_list.append(htos.values[0].split())
                ret_list.append(sub.values[0].split())
            else:
                ret_list.extend(htos.values[0].split())

            return ret_list
        elif hs == None and test_len > 1:
            if not d_con:
                assert sub.str.count(',').values[0] == (
                    htos.str.count(',').values[0]
                    ), "Not equal number of hto and donor names"
                ret_list.append(htos.tolist())
                ret_list.append(sub.tolist())
            else:
                ret_list.extend(htos.tolist())
            
            return ret_list
        elif hs != None and test_len == 1 and \
            htos.str.count(hs).values[0] > 1:
            if not d_con:
                assert sub.str.count(hs).values[0] == (
                    htos.str.count(hs).values[0]
                ), \
                f"""Not equal number of hto and donor names \
                with sep: {hs}"""
                ret_list.append(htos.values[0].split(hs))
                ret_list.append(sub.values[0].split(hs))
            else:
                ret_list.extend(htos.values[0].split(hs))

            return ret_list
        elif hs != None and test_len > 1:
            raise ValueError(
                f"After subsetting sample {s_name} from the wet lab "
                f"file, there are multiple rows ({test_len}) "
                "of HTOs for this sample while a separator value is "
                "also provided."
                )
        elif hs != None and test_len == 1 and htos.str.count(hs).values[0] == 1:
            raise ValueError(
                f"Either the given separator {hs} is wrong or the sample "
                f"{s_name} has incomplete HTO values in the wet lab file "
                )
        else:
            raise ValueError(
                "Something is wrong with the given input(s):\n\twet lab "
                f"file:\n\tsample: {s_name}\n\tHTO-separator: "
                f"{hs}"
                )

    # If 'multiple' HTOs are used for identification
    elif isinstance(hto_col, list):
        # If hto_col is a list of headers
        try:
            return wet_lab_df[hto_col].to_numpy().tolist()
        # If hto_col is a list of integers representing columns
        except:
            return wet_lab_df.iloc[:, [hto_col]].to_numpy().tolist()


    else:
        raise ValueError(
            "Unexpected typing for the columns! Expected either str "
            f"or list but got {type(hto_col)}"
            )


def demux_stats(demux_freq: pd.Series, demux_name: str) -> list[str]:
    r"""
    """

    temp_df=[]
    # Remaining cells after demultiplexing (for each demux method)
    remain_cell = demux_freq.value_counts()
    prop_dict = remain_cell[(remain_cell.index != "Doublet") & \
        (remain_cell.index != "Negative") & \
        (remain_cell.index != "Not Present")]
    od = ord_dict(sorted(prop_dict.items()))
    a = ""
    for k, v in od.items():
        b = str(k) + ": " + str(v) + ", "
        a += b

    temp_df.append(( f'After {demux_name} demultiplexing #cells', 
                    a.strip()[:-1] ))

    return temp_df


# Assumes similar construct like the HTOs
# def parse_subids(wet_lab_df, col_val, fname, s_name, hs=None) -> list:

#     if isinstance(col_val, str):
#         sub = wet_lab_df[col_val]
#         test_len = len(sub)
#         # No command-line params and inference that all HTOs are present 
#         # in one row separated by ","
#         if hs == None and test_len == 1 and sub.str.count(',').values[0] > 1:
#             return sub.values[0].split(',')
#         # No command-line params and inference that all HTOs are present 
#         # in one row separated by whitespaces    
#         elif hs == None and test_len == 1 and len(sub.split()) > 1:
#             return sub.values[0].split(',')
#         elif hs == None and test_len > 1:
#             return sub.tolist()
#         elif hs != None and test_len == 1 and sub.str.count(hs).values[0] > 1:
#             return sub.values[0].split(hs)
#         elif hs != None and test_len > 1:
#             raise ValueError(
#                 f"After subsetting sample {s_name} from the wet lab file "
#                 f", there are multiple rows ({test_len}) of HTOs for "
#                 "this sample while a separator value is also provided."
#             )
#         elif hs != None and test_len == 1 and sub.str.count(hs).values[0] == 1:
#             raise ValueError(
#                 f"Either the given separator {hs} is wrong or the sample "
#                 f"{s_name} has incomplete donor ids in the wet lab file"
#                 )
#         else:
#             raise ValueError(
#                 "Something is wrong with the given input(s):\n\twet lab "
#                 f"file:\n\tsample: {s_name}\n\tHTO-separator: {hs}"
#                 )

def get_donor_info(hto_df: pd.DataFrame, pool_info_df: pd.DataFrame, 
    sep: str, col_list: list):
    r"""Return donor information for each cell barcode for multi-HTO setup.

    This function returns a pandas series containing demultiplexed donor
    info according to the data contained in the wet lab info file. This 
    function is made specially for multi-HTO setup.

    Parameters
    ----------
    hto_df
        A series of cell barcodes from gene count matrix
    pool_info_df
        Subset of wet lab file containing multi-HTO information and \
        SubID (donor IDs)

    col_list
        List of column names in the wet lab file in the sequence: pool name, \
        HTO names (separated by 'sep'), HTO barcode, donor info

    Returns
    -------
    pd.Series
        Contains donor IDs with cell barcodes as index

    """

    # mask=np.column_stack([f for f in ])
    pass



# calico_solo demultiplexing function----------------------------------------
def ret_htos_calico_solo(bcs: pd.Series, df_s: pd.DataFrame, samp: str,
    sep: Optional[str], col_list: list[str, str], dem_cs: pd.Series, 
    donor_convert: bool, hto_count: int, multi_hto_setp: bool
    ) -> list[pd.Series, pd.Series, int, int]:
    r"""Return HTO information and classification for each cell barcode.

    This function returns a 2 pandas series representing donor IDs and 
    HTO name (used for calico_solo) and the number of doublets and 
    negatives identified.

    Parameters
    ----------
    bcs
        A pd series of cell barcodes from gene count matrix
    df_s
        Wet lab file containing HTO information and SubID (donor IDs) \
        for each pool
    samp
        Pool name (present in `df_s` file)
    sep
        Separator used if all HTOs and donors are present in one row \
        or if multi-HTO setup otherwise None
    col_list
        List of column names (first val HTO, second val SubID)
    dem_cs
        A pd series with cell barcodes as index and "HTO classification" \
        (solo output)
    donor_convert
        If donor names have to be converted from the names used in \
        calico_solo (hashsolo) demultiplexing method
    hto_count
        If run for multi-HTO setup this indicates the position of HTO \
        in the sequence
    multi_hto_setp
        True for multi-HTO setup

    Returns
    -------
    pd.Series
        Contains donor IDs with cell barcodes as index
    pd.Series
        Contains HTO name with cell barcodes as index
    int
        number of doublets
    int
        number of negatives.

    """

    # List of htos and don ids from the wet lab spreadsheet (respectively)
    # the indices of these 2 correspond to each other
    if not donor_convert:
        hto_l, subid_l = parse_file(df_s, col_list, samp, sep, donor_convert)
    else:
        hto_l = parse_file(df_s, col_list, samp, sep, donor_convert)

    # List of barcodes
    # barc_l = []
    # SubID from Shan's csv file
    ret_samp = []
    # List of HTOs as HTO1, HTO2, etc
    hash_n = []
    # Doublets' count
    doublet_n = 0
    # Negatives' count
    negative_n = 0

    for bc in bcs:
        if bc in dem_cs.index:
            hto_n = (
                dem_cs[dem_cs.index == bc]
                .values[0]
                )
            
            if hto_n == 'Doublet':
                # barc_l.append(bc)
                hash_n.append(hto_n)
                ret_samp.append(hto_n)
                doublet_n += 1

            elif hto_n == 'Negative':
                # barc_l.append(bc)
                hash_n.append(hto_n)
                ret_samp.append(hto_n)
                negative_n += 1

            else:
                # barc_l.append(bc)
                hash_n.append(hto_n)
                # If the subID doesn't match the HTO value
                if not donor_convert:
                    try:
                        ret_samp.append(subid_l[hto_l.index(hto_n)])
                    except:
                        ret_samp.append(hto_n+"_NA")
                else:
                    ret_samp.append(hto_n)                    

        else:
            # barc_l.append(bc)
            hash_n.append("Not Present")
            ret_samp.append("Not Present")


    ret_samp = pd.Series(ret_samp, index=bcs.index)
    hash_n = pd.Series(hash_n, index=bcs.index)
    # ser_s = pd.DataFrame({'Sample':ret_samp, 'HTO':hash_n}, index=barc_l)

    return [ret_samp, hash_n, doublet_n, negative_n] #zip(barc_l, samp_n)



# def ret_samp_names(y: pd.Series, df_info: pd.DataFrame) -> str:
#     #print(len(y), len(bc))
#     if y in df_info.index:
#         return df_info.loc[df_info.index == y, 'Sample'].values[0]
    
#     else:
#         return "Not Present"



# def ret_hto_number(y: pd.Series, df_info: pd.DataFrame) -> str:
#     #print(len(y), len(bc))
#     if y in df_info.index:
#         return df_info.loc[df_info.index == y, 'HTO'].values[0]

#     else:
#         return "Not Present"


# ---------------------------------------------------------------------------
# Calico_solo execution------------------------------------------------------
# def demux_by_calico_solo(obs_index: pd.Index, df_s: pd.DataFrame):
def demux_by_calico_solo(bcs: pd.Series, df_s: pd.DataFrame, samp: str,
    sep: str, col_list: list[str, str], dem_cs: pd.Series, 
    donor_convert: bool, hto_count: int, multi_hto_sep: str = ""
    ) -> list[pd.Series, pd.Series, list[str]]:
    r"""Main function for classification by calico_solo.

    This function assigns calico_solo classification using another
    function :func:`demultiplex_helper_funcs.ret_htos_calico_solo`.

    """

    # Will contain demultiplexing stats
    temp_df=[]

    # hto_tags_cs = ret_htos_calico_solo(bcs, df_s, samp, sep, col_list, dem_cs)
    SubID_cs, hasht_n_cs, n_doubs, n_negs = ret_htos_calico_solo(bcs, 
                                            df_s, samp, sep, col_list, 
                                            dem_cs, donor_convert)
    
    # Create obs columns in adata to represent the SubID as assigned by 
    # calico solo and its associated hastag number
    # Now DEPRACATED

    # SubID_cs = bcs.to_series().apply(ret_samp_names, args=(hto_tags_cs[0], ))
    # hasht_n_cs = bcs.to_series().apply(ret_hto_number, args=(hto_tags_cs[0], ))

    # Force convert all values to str
    SubID_cs = SubID_cs.apply(str)
    hasht_n_cs = hasht_n_cs.apply(str)
    
    # Save doublets and negatives info from calico solo
    temp_df.append(('Doublets #cells_cs', n_doubs))
    temp_df.append(('Negative #cells_cs', n_negs))

    # Remaining cells after demultiplexing (for each demux method)
    temp_df.extend(demux_stats(SubID_cs, "cs"))

    return [SubID_cs, hasht_n_cs, temp_df]


# ---------------------------------------------------------------------------

# For adding vireoSNP classfication------------------------------------------
def set_don_ids(x: str) -> str:
    r"""Change naming conventions of Vireo

    Use this function to change the naming convetion used in the vireo 
    output - donor_ids.tsv - (generally to make this similar to that 
    of calico_solo/hashsolo but also suits to beautify donor names so 
    as to make it feasible to be classified by using a converter file.

    Parameters
    -----------
    x
        A string representing a donor classification from vireo

    Returns
    -------
    str
        The 'changed' classification
    """

    # Example:
    # One of vireo output - summary.tsv (NOTE: we change it per cell
    # from donor_ids.tsv)
    #      
    # Var1	Freq
    # donor0	82
    # xyzabc	4807
    # qqqqqq	3229
    # nnnnnn	2898
    # donor4	4047
    # donor5	3835
    # doublet	1542
    # unassigned	285
    # 
    # A conversion file that has:
    # vir_out donor_name
    # nnnnnn C1234
    # xyzabc A567
    # qqqqqq A1000
    # 
    # then this function can change the donors in vireo outputs as follows
    # xyzabc -> A567
    # qqqqqq -> A1000
    # nnnnnn -> C1234
    # doublet -> Doublet
    # unassigned -> Negative
    # donor4 -> donor4
    # donor5 -> donor5


    if x == 'doublet':
        return 'Doublet'
    elif x == 'unassigned':
        return 'Negative'
    else:
        return x


# Project-dependent
# def get_don_ids(x: str, t_df: pd.DataFrame) -> str:
#     r""" Coverting the donor names of vireo output

#     Parameters
#     ----------
#     x
#         A string representing a donor
#     t_df
#         A converter file that contains the donor and its converted

#     Returns
#     -------
#     str
#         String representing converted donor name
#     """
#     # Example:
#     # One of vireo output - summary.tsv (NOTE: we change it per cell
#     # from donor_ids.tsv)
#     #      
#     # Var1	Freq
#     # donor0	82
#     # A_xyzabc	4807
#     # A_qqqqqq	3229
#     # A_nnnnnn	2898
#     # donor4	4047
#     # donor5	3835
#     # doublet	1542
#     # unassigned	285
#     # 
#     # A conversion file that has:
#     # vir_out donor_name
#     # nnnnnn C1234
#     # xyzabc A567
#     # qqqqqq A1000
#     # 
#     # then this function will change the donors in vireo output
#     # NOTE: remember to clean with the function 'set_don_ids' 
#     # for this example

#     try:
#         return (
#             t_df.loc[t_df["primary_genotype"] == x, "SubID"].values[0]
#             )
#     except:
#         return x


def ret_subj_ids(ser: list, t_df: pd.DataFrame) -> pd.DataFrame:
    r""" Returns vireo demux stats

    This function returns extra stats from vireo demux output

    Parameters
    ----------
    ser
        A pd series of cell barcodes from gene count matrix
    t_df
        Vireo output (donor_ids.tsv)

    Returns
    -------
    pd.DataFrame
        A dataframe with extra stats
    """
    headers = ["Subj_ID", "prob_max", "prob_doublet"]
    ret_df_l = []
    for x in ser:
        try:
            ret_df_l.append((
                t_df.loc[t_df["cell"].str.strip() == x.strip(), headers]
                .values.flatten().tolist()
                ))
        except:
            ret_df_l.append(["Not Present", "NA", "NA"])

    return pd.DataFrame(ret_df_l, columns=headers, index=ser.index)


# ---------------------------------------------------------------------------

# Vireo execution------------------------------------------------------------
def demux_by_vireo(bcs: pd.Series, vir_out_file: str, 
    conv_df: Optional[pd.DataFrame] = None, donor_col: Optional[str] = None,
    conv_col: Optional[str] = None, pool_col: Optional[str] = None,
    pool_name: Optional[str] = None, 
    ) -> tuple[pd.Series, list[str], Optional[pd.Series]]:
    r"""Main function for classification by vireoSNP.

    This function assigns vireo classification. 

    Paramters
    ---------
    bcs
        A pd series of cell barcodes from gene count matrix.
    vir_out_file
        Path to donor_ids.tsv file produced by vireo.
    conv_df
        A converter file to change donor names in the vireo output, if needed.
    donor_col
        Donor names containing column in the converter file that matches the vireo output.
    conv_col
        Column, in the converter file, containing the converted names.
    pool_col
        Column, in the converter file, containing the pool names.
    pool_name
        Pool name.

    Returns
    -------
    pd.Series
        Classification by vireo per cell
    list
        Demux stats
    pd.Series
        Converted donor names of classification by vireo per cell
    """

    # Will contain demultiplexing stats
    temp_df=[]

    # Read the vireo output
    vir_class = auto_read(vir_out_file)

    vir_class["donor_id"] = vir_class["donor_id"].apply(set_don_ids)
    vir_class["donor_id"] = vir_class["donor_id"].astype('category')
    vir_class.rename(columns={"donor_id":"Subj_ID"}, inplace=True, 
                        errors="raise")
    
    # If conversion of donor names from vireo's output is needed.
    if conv_df is not None:
        if pool_col is not None:
            conv_df = conv_df.loc[ conv_df[pool_col] == pool_name, :]

        map_dict = pd.Series(conv_df[conv_col].values, index=conv_df[donor_col]).to_dict()
			
        # If conversion not found, keep old values
        for v in vir_class["donor_id"].cat.categories:
            if v not in map_dict:
                map_dict[v] = v
        
        vir_class['converted_ID'] = vir_class['Subj_ID'].map(map_dict)

        # DEPRACATED
        # vir_class['Subj_ID'] = (
        #     vir_class['donor_id'].apply(get_don_ids, args=(conv_df,))
        #     )
        # del vir_class['donor_id']


    get_df = ret_subj_ids(bcs, vir_class)

    try:
        n_doubs = vir_class['Subj_ID'].value_counts()['Doublet']
    except:
        n_doubs = 0
    try:
        n_negs = vir_class['Subj_ID'].value_counts()['Negative']
    except:
        n_negs = 0

    # Save doublets and negatives info from calico solo
    temp_df.append(('Doublets #cells_vs', n_doubs))
    temp_df.append(('Negative #cells_vs', n_negs))

    temp_df.extend(demux_stats(get_df["Subj_ID"], "vs"))

    if conv_df is not None:
        return [get_df["Subj_ID"], temp_df, get_df["converted_ID"]]
    else:
        return [get_df["Subj_ID"], temp_df, None]
# ---------------------------------------------------------------------------