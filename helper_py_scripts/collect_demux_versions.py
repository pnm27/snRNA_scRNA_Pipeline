#!/usr/bin/env python3
"""Table of the demultiplex versions used for each pool's final h5ad.

Reads the basic stats written by demux_final.py for each pool and writes
one row per pool: the demultiplex versions (<tool>:<run_tag>), the tools
and run_tags, whether the swap correction was applied, the number of
donors and of cells passing QC.

Example:
python3 collect_demux_versions.py --pools P1 P2 \
    --basic_stats final/P1_basic_stats.tsv final/P2_basic_stats.tsv \
    -o final/final_demux_versions.tsv
"""

import argparse, os
import pandas as pd


# Column of the output -> row of the basic stats
ROWS = {
    'versions': 'Final demultiplex versions',
    'swap_corrected': 'Swap corrected',
    'n_donors': 'Donors in final',
    'n_cells_qc_pass': 'Cells passing mito and basic filter threshold',
}


def get_argument_parser():
    parser = argparse.ArgumentParser(description="Collect the demultiplex "
        "versions of the final h5ads.")
    parser.add_argument('--pools', nargs='+', required=True)
    parser.add_argument('--basic_stats', nargs='+', required=True,
        help="Basic stats of each pool (same order as --pools).")
    parser.add_argument('-o', '--output', required=True)
    return parser


def pool_row(pool: str, fname: str) -> dict:
    stats = pd.read_csv(fname, sep='\t', dtype=str, keep_default_na=False)
    stats = dict(zip(stats.iloc[:, 0], stats.iloc[:, 1]))
    row = {'pool': pool}
    row.update({col: stats.get(name, '') for col, name in ROWS.items()})

    labels = [v.partition(':') for v in row['versions'].split(',') if v]
    row['tools'] = ','.join(t for t, _, _ in labels)
    row['run_tags'] = ','.join(r for _, _, r in labels)
    return row


def main():
    args = get_argument_parser().parse_args()
    if len(args.pools) != len(args.basic_stats):
        raise ValueError("--pools and --basic_stats need the same number "
            "of values")

    df = pd.DataFrame(
        [pool_row(p, f) for p, f in zip(args.pools, args.basic_stats)],
        columns=['pool', 'versions', 'tools', 'run_tags'] + list(ROWS)[1:],
    )
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
