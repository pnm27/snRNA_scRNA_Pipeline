#!/usr/bin/env python3

import pandas as pd, warnings, os, sys, regex
from collections import OrderedDict as ord_dict
from typing import Union, Optional, Any, Literal # Need verion > 3.5


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

# Adding Annotations & validating them --------------------------------------

# EXTRA ANNOTATIONS: TRANSFORMATION
# OLD STYLE
# def apply_pythonic(value, operation, old, new, start, end, by, index):
#     """Apply built-in string methods with optional args."""
#     if operation == "join" and (by is None or by == ' '):
#         method = getattr(' ', operation)
#         return method(value)
#     elif operation == "join" and by is not None:
#         method = getattr(by, operation)
#         return method(value)
#     elif operation == "select_index":
#         if index is not None:
#             return value[index]
#     elif operation == "select_range":
#         if start is None:
#             return value[:end]
#         elif end is None:
#             return value[start:]
#         else:
#             return value[start:end]
#     elif operation == "split":
#         method = getattr(value, operation)
#         return method(' ') if (by is None or by == ' ') else method(by)
#     elif operation == "replace":
#         method = getattr(value, operation)
#         return method(old, new)

from dataclasses import dataclass
# from typing import List, Any


@dataclass
class ColumnResult:
    column_header_h5ad: str
    column_value: str
    column_header_logs: list[str]


from functools import cached_property


class ParsedColumns:
    def __init__(self, results: list[ColumnResult]):
        self.results = results

    @cached_property
    def by_h5ad(self) -> dict[str, str]:
        return {
            r.column_header_h5ad: r.column_value
            for r in self.results
        }

    @cached_property
    def by_logs(self) -> dict[str, str]:
        return {
            tuple(r.column_header_logs): r.column_value
            for r in self.results
        }
    
    

def apply_pythonic(value, t):
    r"""Apply regex substitution.
    
    Paramters
    ---------
    value
        Any value on which python operation(s) needs to be applied
    t
        A dict containing different set of values for different
        python operations

    Returns
    -------
    str
        Output of regex operation
    """

    op = t["operation"]

    if op == "split":
        return value.split(t.get("by", " "))

    elif op == "join":
        return t.get("by", " ").join(value)

    elif op == "select_index":
        return value[t["index"]]

    elif op == "select_range":
        return value[t.get("start"): t.get("end")]

    elif op == "replace":
        return value.replace(t["old"], t["new"])
    

def apply_regex(value: Any, t: dict[str, Any]) -> str:
    r"""Apply regex substitution.
    
    Paramters
    ---------
    value
        Any value on which the regex pattern needs to be applied
    t
        A dict containing 'pattern' and 'replacement' (including flags)

    Returns
    -------
    str
        Output of regex operation
    """
    pattern = t["pattern"]
    replacement = t["replacement"]
    flags = t.get("flags")
    if flags:
        return regex.sub(pattern, replacement, value, flags=flags)
    return regex.sub(pattern, replacement, value)


# EXTRA ANNOTATIONS: PROCESS EACH COLUMN
def process_columns(config: dict[str, Any], 
    pool_name: str, df: pd.DataFrame) -> ParsedColumns:
    r"""Process annotations by reading a JSON recipe.

    This function provides an output, which can be used to 
    annotate sample-level (pool-level) metrics using a template (JSON)
    and parsing a file containing relevant information.

    Paramters
    ---------
    config
        A dict formed from JSON parsing containing recipes for adding 
        annotations.
    pool_name
        Pool name. 
    df
        A pandas dataframe containing relevant annotations.
    Returns
    -------
    ParsedColumns
        An instance containing the results

    See Also
    --------
    :func: ~apply_regex : Regex function can be used, if JSON input needs it.
    :func: ~apply_pythonic : Python function can be used, if JSON input needs it.
    """
    result = []

    for col in config["columns"]:
        current_value = col["source_value"]
        current_value = pool_name if current_value == 'args.p' else current_value
        # Use Wet lab spreadsheet
        if "lookup_column" in col and 'wet_lab.' in col["lookup_column"]:
            col2use = col["value_column"].split('.')[1]
            lookup_column = col["lookup_column"].split('.')[1]
            current_value = df.loc[df[lookup_column] == current_value, col2use].values[0]

        if "transformations" in col:
            for t in col["transformations"]:
                if t["type"] == "pythonic":
                    current_value = apply_pythonic(current_value, t)
                elif t["type"] == "regex":
                    current_value = apply_regex(current_value, t)

        result.append(
            ColumnResult(
                column_header_h5ad=col["name"],
                column_value=current_value,
                column_header_logs=col.get("columnInLogs", [])
            )
        )

    return ParsedColumns(result)


# get demultiplex_paths
def get_demux_paths(config: dict[str, Any]) -> dict[str, str]:
    r""" Function that returns directories with demultiplex results

    This function extracts all the directories containing demultiplexing
    data and returns them.

    Paramters
    ---------
    config
        A dict formed from JSON parsing containing recipes for adding 
        annotations.

    Returns
    -------
    dict
        A dict containing all directories that contain demultiplexed 
        info.
    """
    
    result = {}
    if 'demultiplex_paths' not in config:
        raise warnings.warn("No demultiplex paths provided! "
            "Can't write swap corrected statistics!!!", UserWarning)
    else:
        result.update(config['demultiplex_paths'])

    return result


# PROCESS SWAP CORRECTION
def process_swap_correction(config: dict[str, Any], 
    swap_df: pd.DataFrame, pool_name: str, 
    demux_paths: dict[str, str]) -> str:
    r"""This function returns the location of demultiplexing data

    This function uses a swap correction file to look up specific
    pools and extract the relevant directory containing the correct
    version of demultiplexing run.

    Paramters
    ---------
    config
        A dict formed from JSON parsing containing recipes for adding 
        annotations.
    swap_df
        A pandas dataframe containing swap corrected results.
    pool_name
        Pool name. 
    df
        A pandas dataframe containing relevant annotations.

    Returns
    -------
    str
        A string returning demultiplexing directory.

    See Also
    --------
    :func: ~get_demux_paths : This is utilized to get all possible
    demultiplex output containing directories.
    """

    result = ""
    lookup_columns = [
        'pool_column', 
        'demultiplex_version_column'   
    ]
    
    if 'swap_correction_df' not in config:
        raise warnings.warn("No swap_correction metrics provided! "
            "Can't write swap corrected statistics!!!", UserWarning)
    else:
        if all( f in config['swap_correction_df'] for f in lookup_columns ):
            pool_col = config['swap_correction_df']['pool_column']
            dem_col = config['swap_correction_df']['demultiplex_version_column']
            dem_val = swap_df.loc[swap_df[pool_col] == pool_name, dem_col].values[0]
            result = demux_paths[dem_val]
            
        else:
            raise warnings.warn("Given columns not douns in the swap_correction DF! "
            "Check their names!!!", UserWarning)
    
    return result

# ---------------------------------------------------------------------------

# UPDATED write_logs from update_logs.py

import logging
import pandas as pd

# ---------------------------------------------------------------------------
# File-loading strategies, keyed by sub_prog prefix.
# Each loader receives all_files_dict and returns a DataFrame ready for lookup.
# ---------------------------------------------------------------------------
_LOADERS = {
    "REG": lambda d: pd.read_csv(
        d["STAR_final"],
        names=["cols", "vals"],
        delimiter=r"|",
        skiprows=[7, 22, 27, 34],
    ).assign(
        cols=lambda df: df["cols"].str.strip(),
        vals=lambda df: df["vals"].str.strip(),
    ),
    "GC": lambda d: pd.read_csv(
        d["PICARD_GC"], sep="\t", skiprows=6
    ),
    "RNASEQMETRIC": lambda d: pd.read_csv(
        d["PICARD_RNASeq"], sep="\t", nrows=1, skiprows=6
    ),
    "GENE_FEATURE":     lambda d: get_df(d["Gene_Features"]),
    "GENE_SUMM":        lambda d: pd.read_csv(d["Gene_Summary"],    names=["cols", "vals"]),
    "GENEFULL_FEATURE": lambda d: get_df(d["GeneFull_Features"]),
    "GENEFULL_SUMM":    lambda d: pd.read_csv(d["GeneFull_Summary"], names=["cols", "vals"]),
    "BARCODE_STATS":    lambda d: get_df(d["Barcodes_stats"]),
    "DEMUX":            lambda d: pd.read_csv(
        d["Demultiplex_stats"], names=["cols", "vals"], skiprows=1, sep="\t"
    ),
}

# sub_progs that share the DEMUX loader
_DEMUX_PREFIXES = {"DEMUX", "DEMUX_CS", "DEMUX_VS"}


def _load_file(sub_prog, all_files_dict, cache):
    """
    Return (and cache) the DataFrame for a given sub_prog.
    Raises RuntimeError if the underlying file cannot be read.
    """
    loader_key = "DEMUX" if sub_prog in _DEMUX_PREFIXES else sub_prog
    if loader_key not in cache:
        if loader_key not in _LOADERS:
            raise ValueError(
                f"No file loader defined for sub_prog '{sub_prog}'. "
                "Add an entry to _LOADERS or check your map file."
            )
        try:
            cache[loader_key] = _LOADERS[loader_key](all_files_dict)
        except (FileNotFoundError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            raise RuntimeError(
                f"Could not read source file for sub_prog='{sub_prog}': {e}"
            ) from e
    return cache[loader_key]


def _lookup_value(temp_df, sub_prog, val, mapper):
    """
    Look up a single value from a loaded DataFrame using the mapper.
    Returns the value, or "" if the mapper entry or column is absent.
    Trailing commas are stripped (demux legacy format).
    """
    try:
        mask = (mapper["curr_val"] == val) & (mapper["sub_prog"] == sub_prog)
        mapped_col = mapper.loc[mask, "val_in_log"].values[0]
    except IndexError:
        logging.warning("No mapper entry for sub_prog=%s val=%s — writing empty", sub_prog, val)
        return ""

    try:
        # Key/col-based lookup (GC, RNASEQMETRIC)
        if sub_prog in {"GC", "RNASEQMETRIC"}:
            result = temp_df.loc[0, mapped_col]
        # Row-based lookup: find the row whose 'cols' column matches mapped_col
        else:
            result = temp_df.loc[temp_df["cols"] == mapped_col, "vals"].values[0]
    except (KeyError, IndexError):
        logging.warning(
            "Mapper entry found but value absent in file: sub_prog=%s val=%s mapped_col=%s — writing empty",
            sub_prog, val, mapped_col,
        )
        return ""

    # Normalise: strip trailing commas, clean up floats
    if isinstance(result, str) and result.endswith(","):
        result = result[:-1]
    elif isinstance(result, float) and not result.is_integer():
        result = round(result, 3)

    # REG values use '/' as a separator
    if sub_prog == "REG" and isinstance(result, str):
        result = result.replace(" ", "/")

    return result


def write_logs(big_df, mapper, all_files_dict, no_progs, *args):
    """
    Build one output row for a single sample.

    Parameters
    ----------
    big_df         : DataFrame with MultiIndex columns (prog, sub_prog, val)
    mapper         : DataFrame mapping curr_val -> val_in_log
    all_files_dict : dict of {key: filepath} for all input files
    no_progs       : list of sub_prog prefixes to skip
    *args          : leading annotation values (e.g. sample name)

    Returns
    -------
    list of values in column order
    """
    new_row = list(args)
    file_cache = {}   # Populated lazily; one DataFrame per loader key
    empty_count = 0

    for prog, sub_prog, val in big_df.columns.tolist():

        # LAB columns are annotation placeholders — already handled via *args
        if prog == "LAB":
            continue

        # Skip excluded sub_progs (exact match or prefix match)
        if sub_prog in no_progs or any(sub_prog.startswith(p) for p in no_progs):
            new_row.append("")
            continue

        temp_df = _load_file(sub_prog, all_files_dict, file_cache)
        add_value = _lookup_value(temp_df, sub_prog, val, mapper)

        if add_value == "":
            empty_count += 1

        new_row.append(add_value)

    if empty_count:
        logging.warning(
            "Sample '%s': %d/%d values were empty — check mapper and input files.",
            args[0] if args else "unknown",
            empty_count,
            len(big_df.columns),
        )

    return new_row

# -------------------------------------------------------------------