demultiplex_helper_funcs
========================

.. py:module:: demultiplex_helper_funcs


Attributes
----------

.. autoapisummary::

   demultiplex_helper_funcs._LOADERS
   demultiplex_helper_funcs._DEMUX_PREFIXES


Classes
-------

.. autoapisummary::

   demultiplex_helper_funcs.HashSoloResult
   demultiplex_helper_funcs.ColumnResult
   demultiplex_helper_funcs.ParsedColumns


Functions
---------

.. autoapisummary::

   demultiplex_helper_funcs.auto_read
   demultiplex_helper_funcs.parse_file
   demultiplex_helper_funcs.demux_stats
   demultiplex_helper_funcs.get_donor_info
   demultiplex_helper_funcs._ret_htos_calico_solo
   demultiplex_helper_funcs.demux_by_calico_solo
   demultiplex_helper_funcs.set_don_ids
   demultiplex_helper_funcs.ret_subj_ids
   demultiplex_helper_funcs.demux_by_vireo
   demultiplex_helper_funcs.has_wet_lab_value_column
   demultiplex_helper_funcs.apply_pythonic
   demultiplex_helper_funcs.apply_regex
   demultiplex_helper_funcs.process_columns
   demultiplex_helper_funcs.get_demux_paths
   demultiplex_helper_funcs.process_swap_correction
   demultiplex_helper_funcs.get_df
   demultiplex_helper_funcs._load_file
   demultiplex_helper_funcs._lookup_value
   demultiplex_helper_funcs.write_logs


Module Contents
---------------

.. py:function:: auto_read(fname, lev=1, **kwargs)

.. py:function:: parse_file(wet_lab_df, cols, s_name, hs, d_con)

.. py:function:: demux_stats(demux_freq, demux_name)

   Demux stats


.. py:function:: get_donor_info(hto_df, pool_info_df, sep, col_list)

   Return donor information for each cell barcode for multi-HTO setup.

   This function returns a pandas series containing demultiplexed donor
   info according to the data contained in the wet lab info file. This
   function is made specially for multi-HTO setup.

   :param hto_df:
   :param A series of cell barcodes from gene count matrix:
   :param pool_info_df:
   :param Subset of wet lab file containing multi-HTO information and: SubID (donor IDs)
   :param col_list: List of column names in the wet lab file in the sequence:
                    - pool name
                    - HTO names (separated by 'sep')
                    - HTO barcode
                    - donor info

   :returns: *pd.Series* -- Contains donor IDs with cell barcodes as index


.. py:class:: HashSoloResult

   Container for HashSolo demultiplexing results.


   .. py:attribute:: donor_ids
      :type:  pandas.Series


   .. py:attribute:: hto_names
      :type:  pandas.Series


   .. py:attribute:: n_doublets
      :type:  int


   .. py:attribute:: n_negatives
      :type:  int


.. py:function:: _ret_htos_calico_solo(bcs, df_s, samp, sep, col_list, dem_cs, donor_convert, hto_count, multi_hto_setp)

   Return HTO information and classification for each cell barcode.

   This function returns a 2 pandas series representing donor IDs and
   HTO name (used for calico_solo) and the number of doublets and
   negatives identified.

   :param bcs: A pd series of cell barcodes from gene count matrix
   :param df_s: Wet lab file containing HTO information and SubID (donor IDs) \
                for each pool
   :param samp: Pool name (present in `df_s` file)
   :param sep: Separator used if all HTOs and donors are present in one row \
               or if multi-HTO setup otherwise None
   :param col_list: List of column names (first val HTO, second val SubID)
   :param dem_cs: A pd series with cell barcodes as index and "HTO classification" \
                  (solo output)
   :param donor_convert: If donor names have to be converted from the names used in \
                         calico_solo (hashsolo) demultiplexing method
   :param hto_count: If run for multi-HTO setup this indicates the position of HTO \
                     in the sequence
   :param multi_hto_setp: True for multi-HTO setup

   :returns: *HashSoloResult* -- An instance of the dataclass containing:

             - donor_ids : Sample name indexed by cell barcodes
             - hto_names : HTO name indexed by cell barcodes
             - n_doublets : Total count of doublets
             - n_negatives : Total count of negatives


.. py:function:: demux_by_calico_solo(bcs, df_s, samp, sep, col_list, dem_cs, donor_convert, hto_count, multi_hto_sep = '')

   Assign donor and HTO classifications using HashSolo output.

   This function performs donor assignment using HashSolo (formerly
   calico_solo) classifications generated from HTO demultiplexing.

   :param bcs: Cell barcodes from the gene count matrix.
   :param df_s: Wet lab file containing HTO information and donor IDs
                for each pool.
   :param samp: Pool name present in ``df_s``.
   :param sep: Separator used when all HTOs and donors are present in a
               single row. Set to ``None`` for multi-HTO setup.
   :param col_list: Column names where the first value corresponds to the
                    HTO column and the second value corresponds to the
                    donor ID column.
   :param dem_cs: Series containing HashSolo HTO classifications with
                  cell barcodes as index.
   :param donor_convert: Whether donor names should be converted from the naming
                         convention used by HashSolo.
   :param hto_count: Position of the HTO in the sequence for multi-HTO setup.
   :param multi_hto_sep: Separator used for multi-HTO setups.

   :returns: *tuple* -- A tuple containing:

             ``donor_ids``
                 Donor IDs indexed by cell barcode.

             ``hto_names``
                 HTO names indexed by cell barcode.

             ``demux_stats``
                 Demultiplexing statistics.

   .. seealso::

      :obj:`_ret_htos_calico_solo`
          To generate donor IDs and HTO labels.


.. py:function:: set_don_ids(x)

   Change naming conventions of Vireo

   Use this function to change the naming convetion used in the vireo
   output - donor_ids.tsv - (generally to make this similar to that
   of calico_solo/hashsolo but also suits to beautify donor names so
   as to make it feasible to be classified by using a converter file.

   :param x: A string representing a donor classification from vireo

   :returns: *str* -- The 'changed' classification


.. py:function:: ret_subj_ids(ser, t_df)

   Returns vireo demux stats

   This function returns extra stats from vireo demux output

   :param ser: A pd series of cell barcodes from gene count matrix
   :param t_df: Vireo output (donor_ids.tsv)

   :returns: *pd.DataFrame* -- A dataframe with extra stats


.. py:function:: demux_by_vireo(bcs, vir_out_file, conv_df = None, donor_col = None, conv_col = None, pool_col = None, pool_name = None)

   Main function for classification by vireoSNP.

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

   :returns: * *pd.Series* -- Classification by vireo per cell
             * *list* -- Demux stats
             * *pd.Series* -- Converted donor names of classification by vireo per cell


.. py:function:: has_wet_lab_value_column(config)

   Checks if wet lab file is needed for annotation

   This function evaluates if the Wet Lab file is need for annotations
    as mentioned in the uploaded schema.

   Paramters
   ---------
   config
       A dict formed from JSON parsing containing recipes for adding
       annotations.

   :returns: *bool* -- True if file is needed


.. py:class:: ColumnResult

   .. py:attribute:: column_header_h5ad
      :type:  str


   .. py:attribute:: column_value
      :type:  str


   .. py:attribute:: column_header_logs
      :type:  list[str]


.. py:class:: ParsedColumns(results)

   .. py:attribute:: results


   .. py:property:: by_h5ad
      :type: dict[str, str]



   .. py:property:: by_logs
      :type: dict[str, str]



.. py:function:: apply_pythonic(value, t)

   Apply regex substitution.

   Paramters
   ---------
   value
       Any value on which python operation(s) needs to be applied
   t
       A dict containing different set of values for different
       python operations

   :returns: *str* -- Output of regex operation


.. py:function:: apply_regex(value, t)

   Apply regex substitution.

   Paramters
   ---------
   value
       Any value on which the regex pattern needs to be applied
   t
       A dict containing 'pattern' and 'replacement' (including flags)

   :returns: *str* -- Output of regex operation


.. py:function:: process_columns(config, pool_name, df)

   Process annotations by reading a JSON recipe.

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

   :returns: *ParsedColumns* -- An instance containing the results

   .. seealso::

      :obj:`apply_regex`
          Regex function can be used, if JSON input needs it.
      
      :obj:`apply_pythonic`
          Python function can be used, if JSON input needs it.


.. py:function:: get_demux_paths(config)

   Function that returns directories with demultiplex results

   This function extracts all the directories containing demultiplexing
   data and returns them.

   Paramters
   ---------
   config
       A dict formed from JSON parsing containing recipes for adding
       annotations.

   :returns: *dict* -- A dict containing all directories that contain demultiplexed
             info.


.. py:function:: process_swap_correction(config, swap_df, pool_name, demux_paths)

   This function returns the location of demultiplexing data

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

   :returns: *str* -- A string returning demultiplexing directory.

   .. seealso::

      :obj:`get_demux_paths`
          This is utilized to get all possible demultiplex output containing directories.


.. py:function:: get_df(inp_path)

   This function returns the location of demultiplexing data

   This function uses a swap correction file to look up specific
   pools and extract the relevant directory containing the correct
   version of demultiplexing run.

   Paramters
   ---------
   inp_path
       A dict formed from JSON parsing containing recipes for adding
       annotations.

   :returns: *pd.DataFrame* -- A string returning demultiplexing directory.


.. py:data:: _LOADERS

.. py:data:: _DEMUX_PREFIXES

.. py:function:: _load_file(sub_prog, all_files_dict, cache)

   Return (and cache) the DataFrame for a given sub_prog.
   Raises RuntimeError if the underlying file cannot be read.


.. py:function:: _lookup_value(temp_df, sub_prog, val, mapper)

   Look up a single value from a loaded DataFrame using the mapper.
   Returns the value, or "" if the mapper entry or column is absent.
   Trailing commas are stripped (demux legacy format).


.. py:function:: write_logs(big_df, mapper, all_files_dict, no_progs, *args)

   Build one output row for a single sample.

   :param big_df:
   :type big_df: DataFrame with MultiIndex columns (prog, sub_prog, val)
   :param mapper:
   :type mapper: DataFrame mapping curr_val -> val_in_log
   :param all_files_dict:
   :type all_files_dict: dict of {key: filepath} for all input files
   :param no_progs:
   :type no_progs: list of sub_prog prefixes to skip
   :param \*args:
   :type \*args: leading annotation values (e.g. sample name)

   :returns: *list of values in column order*


