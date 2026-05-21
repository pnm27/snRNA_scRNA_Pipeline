create_inp_cellSNP
==================

.. py:module:: create_inp_cellSNP

.. autoapi-nested-parse::

   Creates a list of effective barcodes for cellSNP run.

   This script takes a count matrix as an input and creates a list of
   'effective' barcodes, dependent on the condition if some cells should
   be removed or not (e.g. if already a demultiplexing tool was run with an
   h5ad column containing classifications and if some of the cells that
   have been classified as a doublet and/or negative cells have to be removed
   then too this script can be used).



Functions
---------

.. autoapisummary::

   create_inp_cellSNP.int_or_none
   create_inp_cellSNP.string_or_none
   create_inp_cellSNP.get_argument_parser
   create_inp_cellSNP.main


Module Contents
---------------

.. py:function:: int_or_none(val)

.. py:function:: string_or_none(val)

.. py:function:: get_argument_parser()

   Generate and return argument parser.


.. py:function:: main()

   Main function


