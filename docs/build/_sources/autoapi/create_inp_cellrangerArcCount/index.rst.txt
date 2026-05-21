create_inp_cellrangerArcCount
=============================

.. py:module:: create_inp_cellrangerArcCount

.. autoapi-nested-parse::

   Creates a metadata file that can be used for cellranger-arc count.

   This script takes as an input all the multiome fastqs (cDNA and ATAC)
   and converts them into a metadata.csv file, which is required for the
   'libraries' parameter in cellranger-arc count command line.



Functions
---------

.. autoapisummary::

   create_inp_cellrangerArcCount.get_argument_parser
   create_inp_cellrangerArcCount.main


Module Contents
---------------

.. py:function:: get_argument_parser()

   Generate and return argument parser.


.. py:function:: main()

   Main function


