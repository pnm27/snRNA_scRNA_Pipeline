create_inp_splitBam_consolidated
================================

.. py:module:: create_inp_splitBam_consolidated


Functions
---------

.. autoapisummary::

   create_inp_splitBam_consolidated.save_df
   create_inp_splitBam_consolidated.get_argument_parser
   create_inp_splitBam_consolidated.main


Module Contents
---------------

.. py:function:: save_df(df, suff, op, overwrite=False)

.. py:function:: get_argument_parser()

   Generate and return argument parser.


.. py:function:: main()

   Main entry function

   This function creates a 2-columned text file with donors as the
   first column with its corresponding barcodes in the second.
   NOT SUPPORTED YET:
   If multiple demultiplexing softwares have annotated the cells
   and one wants to create sep output files for each of them


