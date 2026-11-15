#!/usr/bin/env python3
"""Per-donor bams named by the swap correction file (stand-alone).

The swap correction file has one row per pool and donor: the donor's name
in the split bams ('--donor_column', e.g. the vcf ID) and its name to use
('--name_column', e.g. individualID). Each kept donor's bam is

- linked (hard link; a symlink across filesystems, or '--copy') to
  <out_dir>/<name>.bam, when the name has one bam
- merged (samtools merge) into <out_dir>/<name>.bam, when the name has bams
  of several pools (or several donors of one pool). Each read gets a read
  group of its pool (RG:Z:<pool>): the same cell barcode in two pools is two
  cells

'--keep_apart <column ...>' merges only the bams that also share these
columns (e.g. tissue: <out_dir>/<name>_<tissue>.bam); '--no_combine' never
merges and links each bam to <out_dir>/<name>/<pool>.bam.

Donors with an empty name, or (with '--status_column') a status other than
True / yes / 1, are dropped. Every output and its source bams are listed in
<out_dir>/manifest.tsv; '--dry_run' prints that plan only. Outputs that
already exist are kept unless '--force'.

Example::

    python3 donor_bams.py --swap final_swap_correction.csv \\
        --pool_column pool --donor_column vcf_id --name_column individualID \\
        --status_column final_status \\
        --bams 'split/final/Sample_{pool}-cDNA/{donor}.bam' \\
        --out_dir donor_bams/ --threads 8 --dry_run
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

import pandas as pd


KEPT_STATUS = {'true', 'yes', '1'}
MANIFEST = 'manifest.tsv'
MANIFEST_COLS = ['output', 'action', 'sources', 'pools', 'donors']


def get_argument_parser():
    parser = argparse.ArgumentParser(description="Link or merge per-donor "
        "bams into bams named by the swap correction file.")
    parser.add_argument('--swap', required=True,
        help="Swap correction file (csv, or tsv / txt: tab separated).")
    parser.add_argument('--pool_column', required=True)
    parser.add_argument('--donor_column', required=True,
        help="Donor names of the split bams (e.g. vcf IDs).")
    parser.add_argument('--name_column', required=True,
        help="Names of the output bams (e.g. individualID); empty drops the donor.")
    parser.add_argument('--status_column', default=None,
        help="Donors whose value isn't True (yes, 1) are dropped.")
    parser.add_argument('--bams', required=True,
        help="Path of a split bam, with {pool} and {donor}.")
    parser.add_argument('--out_dir', required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--no_combine', action='store_true',
        help="Never merge: link each bam to <out_dir>/<name>/<pool>.bam.")
    group.add_argument('--keep_apart', nargs='+', default=[],
        metavar='COLUMN', help="Merge only bams that also share these columns.")
    parser.add_argument('--copy', action='store_true',
        help="Copy single bams instead of linking them.")
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--force', action='store_true',
        help="Redo outputs that already exist.")
    parser.add_argument('--dry_run', action='store_true',
        help="Print the plan (manifest) and write nothing.")
    return parser


@dataclass
class Output:
    path: str
    sources: list = field(default_factory=list)   # (pool, donor, bam)
    action: str = ''

    def row(self) -> dict:
        return {
            'output': self.path, 'action': self.action,
            'sources': ','.join(b for _, _, b in self.sources),
            'pools': ','.join(p for p, _, _ in self.sources),
            'donors': ','.join(d for _, d, _ in self.sources),
        }


def read_table(fname: str) -> pd.DataFrame:
    sep = '\t' if fname.endswith(('.tsv', '.txt')) else ','
    return pd.read_csv(fname, sep=sep, dtype=str, keep_default_na=False)


def is_kept(status) -> bool:
    return str(status).strip().lower() in KEPT_STATUS


def safe(name: str) -> str:
    """A name usable as a file name."""
    return re.sub(r'[/\s]+', '_', str(name).strip())


def check_columns(df: pd.DataFrame, args) -> None:
    cols = [args.pool_column, args.donor_column, args.name_column,
            *([args.status_column] if args.status_column else []),
            *args.keep_apart]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        sys.exit(f"{args.swap}: no column(s) {', '.join(missing)} "
                 f"(columns: {', '.join(df.columns)})")
    dup = df.duplicated([args.pool_column, args.donor_column], keep=False)
    if dup.any():
        rows = df.loc[dup, [args.pool_column, args.donor_column]]
        sys.exit(f"{args.swap}: pool and donor repeated in rows:\n"
                 f"{rows.to_string()}")


def plan(df: pd.DataFrame, args) -> tuple[list, list]:
    """Outputs to make and the dropped donors (as Outputs without a path)."""
    names = df[args.name_column].str.strip()
    kept = names != ''
    if args.status_column:
        kept &= df[args.status_column].map(is_kept)

    dropped = [
        Output('', [(r[args.pool_column], r[args.donor_column], '')], 'dropped')
        for _, r in df.loc[~kept].iterrows()
    ]
    df = df.loc[kept].assign(_name=names[kept])

    def source(r):
        return (r[args.pool_column], r[args.donor_column],
                args.bams.format(pool=r[args.pool_column], donor=r[args.donor_column]))

    outputs = []
    if args.no_combine:
        for name, sub in df.groupby('_name', sort=False):
            repeated = sub[args.pool_column].duplicated(keep=False)
            for (_, r), rep in zip(sub.iterrows(), repeated):
                base = (f"{safe(r[args.pool_column])}_{safe(r[args.donor_column])}"
                        if rep else safe(r[args.pool_column]))
                outputs.append(Output(
                    os.path.join(args.out_dir, safe(name), f"{base}.bam"), [source(r)]))
    else:
        keys = ['_name', *args.keep_apart]
        for key, sub in df.groupby(keys, sort=False):
            key = key if isinstance(key, tuple) else (key,)
            base = '_'.join(safe(k) for k in key if str(k).strip() != '')
            outputs.append(Output(
                os.path.join(args.out_dir, f"{base}.bam"),
                [source(r) for _, r in sub.iterrows()]))

    paths = [o.path for o in outputs]
    clash = {p for p in paths if paths.count(p) > 1}
    if clash:
        sys.exit(f"Outputs named the same by different groups: {', '.join(sorted(clash))}")

    link = 'copy' if args.copy else 'link'
    for o in outputs:
        o.action = 'merge' if len(o.sources) > 1 else link
    return outputs, dropped


def run(cmd: list) -> str:
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"Failed: {' '.join(cmd)}\n{res.stderr}")
    return res.stdout


def sq_lines(bam: str) -> list:
    return [l for l in run(['samtools', 'view', '-H', bam]).splitlines()
            if l.startswith('@SQ')]


def link_file(src: str, dst: str, copy: bool) -> None:
    if os.path.lexists(dst):
        os.remove(dst)
    if copy:
        shutil.copy2(src, dst)
        return
    try:
        os.link(src, dst)
    except OSError:     # e.g. another filesystem
        os.symlink(os.path.abspath(src), dst)


def make_single(o: Output, args) -> None:
    bam = o.sources[0][2]
    link_file(bam, o.path, args.copy)
    if os.path.exists(f"{bam}.bai"):
        link_file(f"{bam}.bai", f"{o.path}.bai", args.copy)
    else:
        run(['samtools', 'index', o.path])


def make_merged(o: Output, args) -> None:
    bams = [b for _, _, b in o.sources]
    ref = sq_lines(bams[0])
    for b in bams[1:]:
        if sq_lines(b) != ref:
            sys.exit(f"{b} and {bams[0]} have different references (@SQ); "
                     f"they can't be merged into {o.path}")

    # samtools merge -r names each read group after its input file: link
    # the inputs as <pool>.bam (<pool>_<donor>.bam for donors of one pool)
    tmp = f"{o.path}.tmp_inputs"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    pools = [p for p, _, _ in o.sources]
    inputs = []
    for pool, donor, bam in o.sources:
        rg = f"{safe(pool)}_{safe(donor)}" if pools.count(pool) > 1 else safe(pool)
        inputs.append(os.path.join(tmp, f"{rg}.bam"))
        link_file(bam, inputs[-1], copy=False)

    part = f"{o.path}.part.bam"
    run(['samtools', 'merge', '-f', '-r', '-@', str(args.threads),
         '-o', part, *inputs])
    os.replace(part, o.path)
    shutil.rmtree(tmp)
    run(['samtools', 'index', '-@', str(args.threads), o.path])


def main(argv=None):
    args = get_argument_parser().parse_args(argv)
    if '{pool}' not in args.bams or '{donor}' not in args.bams:
        sys.exit("--bams needs both {pool} and {donor}")

    df = read_table(args.swap)
    check_columns(df, args)
    outputs, dropped = plan(df, args)

    missing = [b for o in outputs for _, _, b in o.sources if not os.path.isfile(b)]
    if missing:
        sys.exit(f"{len(missing)} bam(s) not found:\n" + '\n'.join(missing))

    manifest = pd.DataFrame([o.row() for o in outputs + dropped], columns=MANIFEST_COLS)
    if args.dry_run:
        print(manifest.to_csv(sep='\t', index=False), end='')
        return

    os.makedirs(args.out_dir, exist_ok=True)
    for o in outputs:
        if os.path.exists(o.path) and os.path.exists(f"{o.path}.bai") and not args.force:
            print(f"Exists, kept: {o.path}")
            continue
        os.makedirs(os.path.dirname(o.path) or '.', exist_ok=True)
        print(f"{o.action}: {o.path} <- {', '.join(b for _, _, b in o.sources)}")
        (make_merged if o.action == 'merge' else make_single)(o, args)

    manifest.to_csv(os.path.join(args.out_dir, MANIFEST), sep='\t', index=False)
    print(f"{len(outputs)} bam(s), {len(dropped)} donor(s) dropped; "
          f"see {os.path.join(args.out_dir, MANIFEST)}")


if __name__ == '__main__':
    main()
