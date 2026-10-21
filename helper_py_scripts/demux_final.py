#!/usr/bin/env python3
"""Final h5ad of a pool (rule 'demux_final').

Combines the h5ads written by 'demul_samples.py' (rules 'demux_run_<tool>')
for one pool: the first h5ad is kept as is and the demultiplexing columns of
the others are added to it. 'Donor_final' holds the donor of each cell,
taken from vireo when present (else hashsolo) and renamed by the swap
correction file, if one is given.

Writes two stats files (columns 'Observations' and 'Vals'):
- basic stats: the QC stats (uns['qc_stats'] of the first h5ad) and the
  demultiplex versions used
- demultiplex info: the basic stats and the stats of each tool, as read by
  update_logs.py

Example:
python3 demux_final.py --h5ad solo/wo_gt/P1.h5ad vireo/w_gt/P1.h5ad \
    --stats solo/wo_gt/P1_hashsolo_stats.tsv vireo/w_gt/P1_vireo_stats.tsv \
    --sources solo:wo_gt vireo:w_gt --pool_name P1 -o final/P1_final.h5ad \
    --basic_stats final/P1_basic_stats.tsv --info_out final/P1_final_info.tsv
"""

import argparse, os
import anndata as ad
import pandas as pd
from demultiplex_helper_funcs import auto_read


# Donor columns written by demul_samples.py, in order of preference
DONOR_COLUMNS = ['SubID_vs', 'SubID_vs_cdna', 'SubID_cs']
# 'Donor_final' of the cells whose donor is dropped by the swap correction
REMOVED = "Removed"


def get_argument_parser():
    parser = argparse.ArgumentParser(description="Combine the demultiplexed "
        "h5ads of a pool into its final h5ad.")
    parser.add_argument('--h5ad', nargs='+', required=True,
        help="h5ads of the demultiplexing runs, oldest first.")
    parser.add_argument('--stats', nargs='*', default=[],
        help="Stats files of the demultiplexing tools.")
    parser.add_argument('--sources', nargs='+', required=True,
        help="<tool>:<run_tag> of each h5ad (same order).")
    parser.add_argument('-p', '--pool_name', required=True)
    parser.add_argument('-o', '--output', required=True,
        help="Final h5ad.")
    parser.add_argument('--basic_stats', required=True,
        help="QC stats and the demultiplex versions used.")
    parser.add_argument('--info_out', required=True,
        help="Basic stats and the stats of each tool.")
    parser.add_argument('--donor_column', default=None,
        help="Obs column with the (converted) vireo donor names, if any; "
        "preferred over the columns in DONOR_COLUMNS.")
    parser.add_argument('--swap_correct', default=None,
        help="Swap correction file (csv or tsv).")
    parser.add_argument('--swap_columns', nargs=3, default=None,
        metavar=('pool', 'donor_from', 'donor_to'),
        help="Columns of the swap correction file.")
    return parser


def combine_h5ads(paths: list[str]) -> ad.AnnData:
    """First h5ad, with the obs columns of the others it doesn't have."""
    adata = ad.read_h5ad(paths[0])
    for p in paths[1:]:
        other = ad.read_h5ad(p).obs
        new_cols = [c for c in other.columns if c not in adata.obs.columns]
        adata.obs = adata.obs.join(other[new_cols])
    return adata


def donor_column(adata: ad.AnnData, preferred: str | None) -> str | None:
    for col in [preferred] + DONOR_COLUMNS:
        if col and col in adata.obs.columns:
            return col
    return None


def swap_donor_map(fname: str, cols: list[str], pool: str) -> dict:
    """{donor_from: donor_to} of 'pool' (donor_to is NaN to drop a donor)."""
    pool_col, from_col, to_col = cols
    df = auto_read(fname)
    df = df.loc[df[pool_col].astype(str) == pool]
    return dict(zip(df[from_col].astype(str), df[to_col]))


STAT_COLS = ['Observations', 'Vals']


def qc_stats(adata: ad.AnnData) -> pd.DataFrame:
    if 'qc_stats' not in adata.uns:
        print("No QC stats (uns['qc_stats']) in the first h5ad.")
        return pd.DataFrame(columns=STAT_COLS)
    return pd.DataFrame(adata.uns['qc_stats'])[STAT_COLS]


def main():
    args = get_argument_parser().parse_args()
    if len(args.h5ad) != len(args.sources):
        raise ValueError("--h5ad and --sources need the same number of values")

    adata = combine_h5ads(args.h5ad)
    adata.uns['demux_sources'] = list(args.sources)

    swapped = False
    col = donor_column(adata, args.donor_column)
    if col is not None:
        donors = adata.obs[col].astype(str)
        if args.swap_correct:
            if not args.swap_columns:
                raise ValueError("--swap_columns is needed with --swap_correct")
            donor_map = swap_donor_map(args.swap_correct, args.swap_columns,
                args.pool_name)
            swapped = bool(donor_map)
            donors = donors.map(
                lambda d: d if d not in donor_map
                else (REMOVED if pd.isna(donor_map[d]) else str(donor_map[d]))
            )
        adata.obs['Donor_final'] = pd.Categorical(donors)
    elif args.swap_correct:
        print(f"No donor column in the h5ads of {args.pool_name}; the swap "
            "correction file isn't applied.")

    basic = pd.concat([qc_stats(adata), pd.DataFrame(
        [('Final demultiplex versions', ','.join(args.sources)),
         ('Swap corrected', 'yes' if swapped else 'no')],
        columns=STAT_COLS,
    )])
    info = pd.concat([basic] + [
        pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False)
        for p in args.stats
    ])

    for p in (args.output, args.basic_stats, args.info_out):
        os.makedirs(os.path.dirname(p) or '.', exist_ok=True)
    basic.to_csv(args.basic_stats, sep='\t', index=False)
    info.to_csv(args.info_out, sep='\t', index=False)
    adata.write(args.output)


if __name__ == '__main__':
    main()
