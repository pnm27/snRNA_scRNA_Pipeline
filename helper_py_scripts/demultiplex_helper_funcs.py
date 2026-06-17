#!/usr/bin/env python3

from __future__ import annotations
import pandas as pd, warnings, os, sys
import regex, logging, glob2, logging
from collections import OrderedDict as ord_dict
from typing import Union, Optional, Any # Need verion > 3.10
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from numpy.typing import NDArray


assert sys.version_info >= (3, 10), "This script needs python version >= 3.10!"
# Logger owned by this module
logger = logging.getLogger(__name__)

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
    Demux stats
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
        Pool Info containing file (usually comes from the wet lab)
    Subset of wet lab file containing multi-HTO information and
        SubID (donor IDs)
    col_list
        List of column names in the wet lab file in the sequence: 
        - pool name
        - HTO names (separated by 'sep')
        - HTO barcode
        - donor info

    Returns
    -------
    pd.Series
        Contains donor IDs with cell barcodes as index

    """

    # mask=np.column_stack([f for f in ])
    pass


@dataclass(slots=True)
class HashSoloResult:
    """Container for HashSolo demultiplexing results."""

    donor_ids: pd.Series
    hto_names: pd.Series
    n_doublets: int
    n_negatives: int


# calico_solo demultiplexing function----------------------------------------
def _ret_htos_calico_solo(bcs: pd.Series, df_s: pd.DataFrame, samp: str,
    sep: Optional[str], col_list: list[str, str], dem_cs: pd.Series, 
    donor_convert: bool, 
    hto_count: int, multi_hto_setp: bool
    # ) -> tuple[pd.Series, pd.Series, int, int]:
    ) -> HashSoloResult:
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
    HashSoloResult
        An instance of the dataclass containing:
        
        - donor_ids : Sample name indexed by cell barcodes
        - hto_names : HTO name indexed by cell barcodes
        - n_doublets : Total count of doublets
        - n_negatives : Total count of negatives

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

    # return [ret_samp, hash_n, doublet_n, negative_n] #zip(barc_l, samp_n)
    return HashSoloResult(
        donor_ids=ret_samp,
        hto_names=hash_n,
        n_doublets=doublet_n,
        n_negatives=negative_n,
    )



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
    donor_convert: bool, 
    hto_count: int, multi_hto_sep: str = ""
    ) -> tuple[pd.Series, pd.Series, list[str]]:
    r"""Assign donor and HTO classifications using HashSolo output.

    This function performs donor assignment using HashSolo (formerly
    calico_solo) classifications generated from HTO demultiplexing.

    Parameters
    ----------
    bcs
        Cell barcodes from the gene count matrix.
    df_s
        Wet lab file containing HTO information and donor IDs
        for each pool.
    samp
        Pool name present in ``df_s``.
    sep
        Separator used when all HTOs and donors are present in a
        single row. Set to ``None`` for multi-HTO setup.
    col_list
        Column names where the first value corresponds to the
        HTO column and the second value corresponds to the
        donor ID column.
    dem_cs
        Series containing HashSolo HTO classifications with
        cell barcodes as index.
    donor_convert 
        Whether donor names should be converted from the naming
        convention used by HashSolo.
    hto_count
        Position of the HTO in the sequence for multi-HTO setup.
    multi_hto_sep
        Separator used for multi-HTO setups.

    Returns
    -------
    Returns a tuple with following values in the order:

    donor_ids
        Donor IDs indexed by cell barcode.

    hto_names
        HTO names indexed by cell barcode.

    demux_stats
        Demultiplexing statistics.

    See Also
    --------
    _ret_htos_calico_solo 
        To generate donor IDs and HTO labels.

    """

    # Will contain demultiplexing stats
    temp_df=[]

    # hto_tags_cs = ret_htos_calico_solo(bcs, df_s, samp, sep, col_list, dem_cs)
    results = _ret_htos_calico_solo(bcs, df_s, samp, sep, 
            col_list, dem_cs, donor_convert
            )
    
    
    # Create obs columns in adata to represent the SubID as assigned by 
    # calico solo and its associated hastag number
    # Now DEPRACATED

    # SubID_cs = bcs.to_series().apply(ret_samp_names, args=(hto_tags_cs[0], ))
    # hasht_n_cs = bcs.to_series().apply(ret_hto_number, args=(hto_tags_cs[0], ))

    # Force convert all values to str
    SubID_cs = results.donor_ids.apply(str)
    hasht_n_cs = results.hto_names.apply(str)
    n_doubs = results.n_doublets
    n_negs = results.n_negatives
    
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

def has_wet_lab_value_column(config: dict) -> bool:
    r"""Checks if wet lab file is needed for annotation

    This function evaluates if the Wet Lab file is need for annotations
     as mentioned in the uploaded schema.

    Paramters
    ---------
    config
        A dict formed from JSON parsing containing recipes for adding 
        annotations.
    
    Returns
    -------
    bool
        True if file is needed

    """
    return any(
        col.get("value_column", "").startswith("wet_lab.")
        for col in config.get("columns", [])
        if "value_column" in col
    )


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


@dataclass
class ColumnResult:
    column_header_h5ad: str
    column_value: str
    column_header_logs: list[str]


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
    pool_name: str, df: pd.DataFrame | None) -> ParsedColumns:
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
    apply_regex : Regex function can be used, if JSON input needs it.
    apply_pythonic : Python function can be used, if JSON input needs it.

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


# ----------------------------------------------------------------------------
# For creating All_logs file containing all metrics --------------------------
def get_filename(loc_dir: str | None, 
    file_struct: str | None, 
    fn: str | None, 
    suffix: str | None) -> list[str]:
    """
    Find matching file(s) in one or more directories.

    Parameters
    ----------
    loc_dir
        Directory path or comma-separated list of directory paths.
    file_struct
        File structure/prefix pattern. If it ends with ``'/'``, it is treated
        as a subdirectory. If empty, the search is based only on ``fn`` and
        ``suffix``.
    fn
        Filename identifier used in the search pattern.
    suffix
        File suffix to match. If ``None``, an empty string is returned.

    Returns
    -------
    str
        Comma-separated list of matching file paths. Missing matches are
        represented by empty entries. Returns an empty string if no valid
        matches are found.
    """

    if suffix is None or not loc_dir:
        return []

    file_struct = file_struct or ""
    fn = fn or ""

    directories = (
        [loc_dir]
        if "," not in loc_dir
        else [d.strip() for d in loc_dir.split(",")]
    )

    results = []

    for directory in directories:
        match = ""

        if file_struct.endswith("/"):
            pattern = Path(directory) / file_struct / f"{fn}*{suffix}"
            matches = glob2.glob(str(pattern))

        elif file_struct == "":
            pattern = Path(directory) / f"{fn}*{suffix}"
            matches = glob2.glob(str(pattern))

        else:
            patterns = [
                Path(directory) / f"{file_struct}*{suffix}",
                Path(directory) / f"{file_struct}*{fn}*{suffix}",
            ]

            matches = []
            for pattern in patterns:
                matches = glob2.glob(str(pattern))
                if matches:
                    break

        if matches:
            match = matches[0]

        results.append(match)


    return results


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
        warnings.warn("""
            No demultiplex paths provided!
            Can't write swap corrected statistics!!!
            """, UserWarning
        )
    else:
        result.update(config['demultiplex_paths'])

    return result


# PROCESS SWAP CORRECTION
def process_swap_correction(config: dict[str, Any], 
    swap_df: pd.DataFrame | None, pool_name: str, 
    demux_paths: dict[str, str],
    logger: logging.Logger | None = None) -> str | None:
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
    demux_paths
        A dict containing all directories that contain demultiplexed 
        info.
    logger
        Logger to use. Defaults to this module's logger.

    Returns
    -------
    str
        A string returning demultiplexing directory.
    None
        When `pool_name` is not found in the `swap_df` or 
        when the required columns are not found in the `swap_df`.

    See Also
    --------
    get_demux_paths : This is utilized to get all possible 
        demultiplex output containing directories.

    """

    result = ""
    lookup_columns = [
        'pool_column', 
        'demultiplex_version_column'   
    ]
    colsFound_msg = (
        f"{pool_name} not present in the swap_correction file! "
        "Skipping writing the demultiplexing info!!!"
    )
    notFound_msg = (
        "No swap_correction metrics provided! "
        "Can't write swap corrected statistics!!!"
    )
    colsNotFound_msg = (
        "Given columns not present in the swap_correction file! " 
        "Check their names!!!"
    )
    if 'swap_correction_df' not in config:
        warnings.warn(notFound_msg, UserWarning)
    else:
        if all( f in config['swap_correction_df'] for f in lookup_columns ):
            pool_col = config['swap_correction_df']['pool_column']
            dem_col = config['swap_correction_df']['demultiplex_version_column']
            try:
                dem_val = swap_df.loc[swap_df[pool_col] == pool_name, dem_col].values[0]
            except IndexError:
                warnings.warn(colsFound_msg, UserWarning)
                logger.warning(colsFound_msg)
                return None
            result = demux_paths[dem_val]
            
        else:
            warnings.warn(colsNotFound_msg, UserWarning)
    
    return result

# ---------------------------------------------------------------------------

# This function is to read files with extension '.stats', which are formatted weirdly
# and return a pandas Dataframe for easy use
def get_df(inp_path) -> pd.DataFrame:
    r"""This function returns the location of demultiplexing data

    This function uses a swap correction file to look up specific
    pools and extract the relevant directory containing the correct
    version of demultiplexing run.

    Paramters
    ---------
    inp_path
        A dict formed from JSON parsing containing recipes for adding 
        annotations.

    Returns
    -------
    pd.DataFrame
        A string returning demultiplexing directory.
    """
    col1 = []
    col2 = []
    with open(inp_path) as f1:
        for line in f1:
            col1.append(line.strip().split()[0])
            col2.append(line.strip().split()[1])
   
    n_df = pd.DataFrame({'cols':col1,'vals':col2})
    return n_df


# UPDATED write_logs from update_logs.py
# ---------------------------------------------------------------------------
# File-loading strategies, keyed by sub_prog prefix.
# Each loader receives all_files_dict and returns a DataFrame ready for lookup.
# ---------------------------------------------------------------------------
_LOADERS = {
    "REG": lambda d: pd.read_csv(
        d["STAR_final"][0],
        names=["cols", "vals"],
        delimiter=r"|",
        skiprows=[7, 22, 27, 34],
    ).assign(
        cols=lambda df: df["cols"].str.strip(),
        vals=lambda df: df["vals"].str.strip(),
    ),
    "GC": lambda d: pd.read_csv(
        d["PICARD_GC"][0], sep="\t", skiprows=6
    ),
    "RNASEQMETRIC": lambda d: pd.read_csv(
        d["PICARD_RNASeq"][0], sep="\t", nrows=1, skiprows=6
    ),
    "GENE_FEATURE":     lambda d: get_df(d["Gene_Features"][0]),
    "GENE_SUMM":        lambda d: pd.read_csv(d["Gene_Summary"][0],    names=["cols", "vals"]),
    "GENEFULL_FEATURE": lambda d: get_df(d["GeneFull_Features"][0]),
    "GENEFULL_SUMM":    lambda d: pd.read_csv(d["GeneFull_Summary"][0], names=["cols", "vals"]),
    "BARCODE_STATS":    lambda d: get_df(d["Barcodes_stats"][0]),
    # FUTURE
    # "DEMUX_CS":            lambda d: pd.read_csv(
    #     d["Demultiplex_stats"], names=["cols", "vals"], skiprows=1, sep="\t"
    # ),
    # "DEMUX_VS":            lambda d: pd.read_csv(
    #     d["Demultiplex_stats"], names=["cols", "vals"], skiprows=1, sep="\t"
    # ),
    "DEMUX":            lambda d: pd.read_csv(
        d["Demultiplex_stats"][0], names=["cols", "vals"], skiprows=1, sep="\t"
    ),
}

# sub_progs that share the DEMUX loader
_DEMUX_PREFIXES = {"DEMUX", "DEMUX_CS", "DEMUX_VS"}
def _load_file(sub_prog: str, 
    all_files_dict: dict[str, list[str]], 
    cache: dict, 
    logger: logging.Logger| None = None) -> pd.DataFrame:
    """
    Return (and cache) the DataFrame for a given sub_prog.

    Parameters
    ----------
    sub_prog
        Sub-program name.
    all_files_dict
        Dictionary containing source file paths.
    cache
        Cache of already-loaded DataFrames.
    logger
        Logger to use. Defaults to this module's logger.

    Returns
    -------
    pandas.DataFrame
        Loaded DataFrame.

    Raises
    ------
    ValueError
        If no loader is defined for ``sub_prog``.
    RuntimeError
        If the source file cannot be read.
    """
    log = logger or logging.getLogger(__name__)

    loader_key = "DEMUX" if sub_prog in _DEMUX_PREFIXES else sub_prog

    if loader_key in cache:
        log.debug(
            "Using cached DataFrame for sub_prog='%s' (loader='%s')",
            sub_prog,
            loader_key,
        )
        return cache[loader_key]

    log.debug(
        "Cache miss for sub_prog='%s'; using loader '%s'",
        sub_prog,
        loader_key,
    )

    if loader_key not in _LOADERS:
        msg = (
            f"No file loader defined for sub_prog '{sub_prog}'. "
            "Add an entry to _LOADERS or check your map file."
        )
        log.error(msg)
        raise ValueError(msg)

    try:
        log.debug("Loading source data for sub_prog='%s'", sub_prog)

        cache[loader_key] = _LOADERS[loader_key](all_files_dict)

        log.info(
            "Successfully loaded source data for sub_prog='%s'",
            sub_prog,
        )

    except (
        FileNotFoundError,
        pd.errors.ParserError,
        pd.errors.EmptyDataError,
    ) as e:

        log.exception(
            "Failed to load source data for sub_prog='%s'",
            sub_prog,
        )

        raise RuntimeError(
            f"Could not read source file for sub_prog='{sub_prog}': {e}"
        ) from e

    return cache[loader_key]


_SKIP_METRICS =[
    "N_DOUBLET_CELLS",
    "N_NEGATIVE_CELLS",
    "N_CELLS_AFTER_DEMUX",
]
def _lookup_value(temp_df: pd.DataFrame, 
    sub_prog: str, 
    metric: str, 
    mapper: pd.DataFrame, 
    logger=None) -> str | int | float:
    """
    Look up a single value from a loaded DataFrame using the mapper.
    Returns the value, or "" if the mapper entry or column is absent.
    Trailing commas are stripped (demux legacy format).

    Parameters
    ----------
    temp_df
        File containing the necessary metrics
    sub_prog
        Sub-program name.
    metric
        Metric name to lookup in `mapper`.
    mapper
        Mapping curr_val -> val_in_log
    logger
        Logger to use. Defaults to this module's logger.

    Returns
    -------
    Returns the looked-up value, or an empty string if not found.

    """
    log = logger or logging.getLogger(__name__)
    try:
        mask = (mapper["curr_val"] == metric) & (mapper["sub_prog"] == sub_prog)
        mapped_col = mapper.loc[mask, "val_in_log"].values[0]
    except IndexError:
        log.warning("No mapper entry for sub_prog=%s val=%s — writing empty", sub_prog, metric)
        return ""

    try:
        # Key/col-based lookup (GC, RNASEQMETRIC)
        if sub_prog in {"GC", "RNASEQMETRIC"}:
            result = temp_df.loc[0, mapped_col]
        # Row-based lookup: find the row whose 'cols' column matches mapped_col
        else:
            result = temp_df.loc[temp_df["cols"] == mapped_col, "vals"].values[0]
    except (KeyError, IndexError):
        if not any(m in metric for m in _SKIP_METRICS):
            log.warning(
                (
                "Mapper entry found but value absent in file: "
                "sub_prog=%s metric=%s mapped_col=%s — writing empty"
                ),
                sub_prog, metric, mapped_col,
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


def write_logs(columns: NDArray, 
    mapper: pd.DataFrame, 
    all_files_dict: dict[str, list[str]], 
    no_progs: list, sample: str, logger=None, 
    processed_data=None):
    """
    Build one output row for a single sample.

    Parameters
    ----------
    columns
        Column names present as prog, sub_prog, val
    mapper
        DataFrame mapping curr_val -> val_in_log
    all_files_dict
        all input files present as key: list of file paths
    no_progs
        list of sub_prog prefixes to skip
    sample
        Sample name (for logging)
    logger
        Logger instance for logging warnings (optional)
    processed_data
        Dict of extra annotations

    Returns
    -------
    list of values in column order
    """
    log = logger or logging.getLogger(__name__)
    # new_row = list(args)
    new_row = []
    file_cache = {}   # Populated lazily; one DataFrame per loader key
    empty_count = 0

    for prog, sub_prog, metric in columns.tolist():

        # LAB columns are always annotated here as they don't require file lookups
        if prog == "LAB":
            new_row.append(processed_data.get((prog, sub_prog, metric), ""))
            continue

        # Skip excluded sub_progs (exact match or prefix match)
        if sub_prog in no_progs or any(sub_prog.startswith(p) for p in no_progs):
            new_row.append("")
            continue

        temp_df = _load_file(sub_prog, all_files_dict, file_cache, logger=log)
        add_value = _lookup_value(temp_df, sub_prog, metric, mapper, logger=log)

        if add_value == "":
            empty_count += 1

        new_row.append(add_value)

    if empty_count:
        log.warning(
            "Sample '%s': %d/%d values were empty — check mapper and input files.",
            sample,
            empty_count,
            len(columns),
        )

    print(f"Finished processing {sample}")

    return new_row

# -------------------------------------------------------------------