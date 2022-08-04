#!/usr/bin/env python

from itertools import repeat
import pandas as pd, numpy as np
import sys, os, re, glob2, argparse
from collections import Counter


#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Demultiplex sample based on hahsolo produced output")

parser.add_argument('matrix_file', help="Path to matrix.mtx.gz")
parser.add_argument('count_matrix', help="Path to store the final count matrix(h5ad)")
parser.add_argument('demux_info', help="Path to store demultiplexing info (tab-separated txt file)")
parser.add_argument('gene_info_file', help="Path to file that contains gene names and ids for annotation (tab-separated txt file)")
parser.add_argument('wet_lab_file', help="Path to file that contains HTO info for each set (either csv or tsv file)")


# Optional parameters
parser.add_argument('--calico_solo', dest='hashsolo_out', help="Path to cached output of hashsolo(h5ad)")
parser.add_argument('--vireo_out', help="Path to cached output of hashsolo(h5ad)")
parser.add_argument('-m', '--max_mito', type=int, help="Max mitochondrial genes(in percent) per cell. Default: 5", default=5)
parser.add_argument('-g', '--min_genes', type=int, help="Min #genes per cell. Default: 1000", default=1000)
parser.add_argument('-c', '--min_cells', type=int, help="Min #cells expressing a gene for it to pass the filter. Default: 10", default=10)
parser.add_argument('--columns', nargs=4, help="List of column names RESPECTIVELY to sample_ID (as present in the spreadsheets - identifies pooled samples), \
    HTO numbers, it's associated barcodes, and Donors/SubIDs (contains each multiplexed donors).", metavar=('sample_ID', 'HTO_name', 'HTO_barcode', 'Sub_ID'),
     default=['unique_sample_ID', 'hashtag', 'ab_barcode', 'SubID'])
parser.add_argument('-s', '--sample_name', help="Name of the sample. Check whether this name coincides with that in this script as well as the one in the wet_lab_file")
parser.add_argument('--hto_sep', help="If, per each sample in the wet lab file (6th positional argument to this script), HTOs are all present in one row separated by some SEP then specify it here. Default: ' '", default=' ')


args = parser.parse_args()

