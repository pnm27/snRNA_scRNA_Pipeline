#!/usr/bin/env python3
"""Donor labels of cells from outside the pipeline, as one table.

Input ('modules: split_bams: split_by: external'): a pooled h5ad whose obs
has the donor and the pool of origin of each cell, or a table (csv, tsv or
txt, optionally gzipped) with these columns. The cell barcodes are made
unique in such files (e.g. 'P12_L1_AAACCTGAGC-1'), so they are transformed
back into the bam's CB tag with 'pythonic'/'regex' steps (see
transform_steps.py); pool names can be transformed the same way, into the
'{pool}' names of the bams.

Output: a tab-separated file with the columns pool, Subj_ID and barcodes,
one row per cell, without doublets, negatives and unassigned cells. The
checkpoint 'create_inp_splitBams' takes the rows of each pool from it.

Usage::

    python3 normalize_external_labels.py pooled.h5ad labels.tsv \\
        --donor_column individualID --pool_column pool_of_origin \\
        --barcode_steps '[{"type": "pythonic", "operation": "split", "by": "_"},
                          {"type": "pythonic", "operation": "select_index", "index": -1}]'
"""

from __future__ import annotations
import argparse, json, os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transform_steps import apply_steps  # noqa: E402

# Labels that aren't donors: demul_samples.py (both tools), the raw
# outputs of the tools and demux_final.py
NOT_DONORS = {'Doublet', 'Negative', 'Not Present', 'doublet', 'unassigned',
              'Removed'}

OUT_COLUMNS = ['pool', 'Subj_ID', 'barcodes']


def read_labels(path: str, donor_column: str, pool_column: str,
                barcode_column: str | None) -> pd.DataFrame:
    """The barcode, donor and pool of each cell in 'path' (h5ad or table)."""
    if path.endswith('.h5ad'):
        import anndata as ad
        obs = ad.read_h5ad(path, backed='r').obs  # obs only, not the matrix
        _require(obs.columns, [donor_column, pool_column], path, 'obs column')
        barcodes = (obs.index if barcode_column is None
                    else obs[_require(obs.columns, [barcode_column], path,
                                      'obs column')[0]])
        df = pd.DataFrame({
            'barcodes': pd.Series(barcodes, dtype=str).to_numpy(),
            'Subj_ID': obs[donor_column].to_numpy(),
            'pool': obs[pool_column].to_numpy(),
        })
    else:
        if barcode_column is None:
            raise ValueError(
                f"{path}: tables need 'barcode_column' (the column with the "
                "cell barcodes); only h5ads can use their obs names")
        name = path[:-3] if path.endswith('.gz') else path
        sep = ',' if name.endswith('.csv') else '\t'
        cols = [barcode_column, donor_column, pool_column]
        header = pd.read_csv(path, sep=sep, nrows=0).columns
        _require(header, cols, path, 'column')
        df = pd.read_csv(path, sep=sep, usecols=cols, dtype=str)
        df = df.rename(columns={barcode_column: 'barcodes',
                                donor_column: 'Subj_ID', pool_column: 'pool'})
    return df[OUT_COLUMNS]


def _require(columns, wanted: list, path: str, what: str) -> list:
    missing = [c for c in wanted if c not in columns]
    if missing:
        raise ValueError(
            f"{path}: no {what}(s) {', '.join(map(repr, missing))}; found: "
            f"{', '.join(map(str, columns))}")
    return wanted


def transform(values: pd.Series, steps: list | None) -> pd.Series:
    """'steps' applied to each value (once per distinct value)."""
    if not steps:
        return values.astype(str)
    uniq = pd.unique(values.astype(str))
    mapping = {v: str(apply_steps(v, steps)) for v in uniq}
    return values.astype(str).map(mapping)


def normalize(df: pd.DataFrame, barcode_steps: list | None = None,
              pool_steps: list | None = None) -> pd.DataFrame:
    """Labels with the bam's barcodes and pool names, donors only."""
    df = df.dropna(subset=['Subj_ID', 'pool'])
    df = df[~df['Subj_ID'].astype(str).isin(NOT_DONORS)].copy()
    df['Subj_ID'] = df['Subj_ID'].astype(str)
    df['barcodes'] = transform(df['barcodes'], barcode_steps)
    df['pool'] = transform(df['pool'], pool_steps)

    dup = df.duplicated(['pool', 'barcodes'], keep=False)
    if dup.any():
        ex = df.loc[dup].head(3).to_string(index=False)
        raise ValueError(
            f"{int(dup.sum())} cells share a barcode within a pool after "
            "the barcode transformations, so their reads can't be told "
            "apart. Check 'barcode_transformations'. E.g.:\n" + ex)
    return df[OUT_COLUMNS].sort_values(['pool', 'Subj_ID'], kind='stable')


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="Donor labels of cells (pooled h5ad or table) as a "
        "pool, Subj_ID, barcodes table, with the bam's barcodes.")
    parser.add_argument('inp', help="Pooled h5ad, or csv/tsv/txt(.gz) table")
    parser.add_argument('output', help="Output table (tab-separated)")
    parser.add_argument('--donor_column', default='individualID',
                        help="Column with the donors (default: individualID)")
    parser.add_argument('--pool_column', required=True,
                        help="Column with the pool of origin of each cell")
    parser.add_argument('--barcode_column', default=None,
                        help="Column with the barcodes (h5ad default: obs names)")
    parser.add_argument('--barcode_steps', type=json.loads, default=None,
                        help="JSON list of steps for the barcodes")
    parser.add_argument('--pool_steps', type=json.loads, default=None,
                        help="JSON list of steps for the pool names")
    return parser


def main():
    args = get_argument_parser().parse_args()
    df = read_labels(args.inp, args.donor_column, args.pool_column,
                     args.barcode_column)
    n_in = len(df)
    df = normalize(df, args.barcode_steps, args.pool_steps)
    if df.empty:
        raise ValueError(f"{args.inp}: no cells with a donor")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, sep='\t', index=False)

    print(f"{len(df)} of {n_in} cells have a donor; per pool:")
    print(df.groupby('pool').agg(cells=('barcodes', 'size'),
                                 donors=('Subj_ID', 'nunique')).to_string())
    print("Example barcodes:", ', '.join(df['barcodes'].head(3)))


if __name__ == '__main__':
    main()
