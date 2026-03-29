#!/usr/bin/env python3

"""Creates per donor h5ads that can be used for downstream processing.

This script takes all h5ad files in a directory, assigns final (swap/gt)
corrected donor names, adds annotations and creates per-donor h5ads.
"""

import anndata as ad, string
import scanpy as sc, pandas as pd, numpy as np
import os, re, sys, argparse
from collections import Counter
from collections import defaultdict
# using datetime module
from demultiplex_helper_funcs import auto_read


# Multiprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def process_file(f, op_png_dir):