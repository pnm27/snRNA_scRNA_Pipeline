#!/usr/bin/env python3
"""Checks of the WASP stage's inputs (rules/wasp.smk), before the long jobs.

samples
    each genotype sample ID of the samples table is in its set's
    haps.h5 (else find_intersecting_snps.py fails hours in)
contigs
    a bam and a set's snp_tab.h5 name chromosomes the same way
    (chr1 vs 1): with no chromosome in common, WASP would find
    no SNP and keep every read
set-samples
    the genotype sample IDs of a set (for bcftools view -S)

The h5 files are read with PyTables (or h5py), as in the WASP module's
python. snp2h5 writes a table 'samples_<chrom>' (column 'name') per
chromosome in haps.h5 and a node per chromosome in snp_tab.h5.

Usage::

    python wasp_check_inputs.py samples --table samples.tsv \\
        --set main=/geno/haps.h5,/geno/snp_tab.h5 -o checked.txt
    python wasp_check_inputs.py contigs --contigs contigs.txt \\
        --snp_tab /geno/snp_tab.h5 --bam donor.bam
    python wasp_check_inputs.py set-samples --table samples.tsv --set main --sets main
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lib.wasp import read_samples_table, set_samples  # noqa: E402


def _decode(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def h5_nodes(path) -> dict:
    """{top-level node name: node data reader} of an h5 file."""
    try:
        import tables
        with tables.open_file(path, 'r') as f:
            out = {}
            for node in f.root:
                name = node._v_name
                if name.startswith('samples_'):
                    out[name] = [_decode(x) for x in node.col('name')]
                else:
                    out[name] = None
            return out
    except ImportError:
        import h5py
        with h5py.File(path, 'r') as f:
            return {name: ([_decode(x) for x in f[name]['name']]
                           if name.startswith('samples_') else None)
                    for name in f.keys()}


def h5_samples(haps) -> set:
    """Genotype sample IDs in a haps.h5 (empty if none could be read)."""
    return {s for name, v in h5_nodes(haps).items()
            if name.startswith('samples_') for s in v or []}


def h5_chromosomes(snp_tab) -> set:
    return {name for name in h5_nodes(snp_tab) if not name.startswith('samples_')}


def missing_samples(table: dict, gset: str, known: set) -> list:
    """'<pool> <donor>: <sample>' of the set's samples not in 'known'."""
    return [f"{pool} {donor}: {s}"
            for (pool, donor), row in sorted(table.items())
            if row['genotype_set'] == gset
            for s in row['samples'] if s not in known]


def common_contigs(bam_contigs, h5_chroms) -> set:
    return set(bam_contigs) & set(h5_chroms)


def check_samples(table_path, sets: dict, need_bam: bool, out=None):
    table = read_samples_table(table_path, need_bam=need_bam, sets=list(sets))
    errors, lines = [], []
    for gset, (haps, snp_tab) in sets.items():
        known = h5_samples(haps)
        n = len(set_samples(table, gset))
        if not known:
            print(f"Warning: no sample names read from {haps} (set {gset}); "
                  "the samples table isn't checked against it", file=sys.stderr)
        else:
            missing = missing_samples(table, gset, known)
            if missing:
                errors.append(
                    f"Genotype set '{gset}' ({haps}) has no sample:\n    "
                    + "\n    ".join(missing))
        chroms = sorted(h5_chromosomes(snp_tab))
        if not chroms:
            errors.append(f"Genotype set '{gset}': no chromosome in {snp_tab}")
        lines.append(f"{gset}\t{n} samples\t{len(chroms)} chromosomes: {','.join(chroms)}")
    if errors:
        raise SystemExit(
            f"WASP samples table {table_path}: the 'sample' column must hold "
            "the genotype IDs of each set (find_intersecting_snps.py's "
            "--samples)\n  " + "\n  ".join(errors))
    if out:
        with open(out, 'w') as f:
            f.write("set\tsamples\tchromosomes\n" + "\n".join(lines) + "\n")


def check_contigs(contigs_path, snp_tab, bam):
    with open(contigs_path) as f:
        contigs = [l.strip() for l in f if l.strip()]
    chroms = h5_chromosomes(snp_tab)
    common = common_contigs(contigs, chroms)
    if not common:
        raise SystemExit(
            f"{bam} and {snp_tab} have no chromosome in common (bam: "
            f"{', '.join(contigs[:5])}...; genotypes: {', '.join(sorted(chroms)[:5])}...). "
            "Do both use the 'chr' prefix (or neither)?")
    return common


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('samples')
    s.add_argument('--table', required=True)
    s.add_argument('--set', action='append', required=True,
                   help="<name>=<haps.h5>,<snp_tab.h5>, once per genotype set")
    s.add_argument('--need_bam', action='store_true')
    s.add_argument('-o', '--output')

    c = sub.add_parser('contigs')
    c.add_argument('--contigs', required=True, help="the bam's chromosomes, one per line")
    c.add_argument('--snp_tab', required=True)
    c.add_argument('--bam', default='the bam')

    g = sub.add_parser('set-samples')
    g.add_argument('--table', required=True)
    g.add_argument('--set', required=True)
    g.add_argument('--sets', nargs='+', required=True, help="all genotype set names")
    g.add_argument('--need_bam', action='store_true')
    args = p.parse_args(argv)

    if args.cmd == 'samples':
        sets = {}
        for item in args.set:
            name, files = item.split('=', 1)
            haps, snp_tab = files.split(',', 1)
            sets[name] = (haps, snp_tab)
        check_samples(args.table, sets, args.need_bam, args.output)
    elif args.cmd == 'contigs':
        check_contigs(args.contigs, args.snp_tab, args.bam)
    else:
        table = read_samples_table(args.table, need_bam=args.need_bam, sets=args.sets)
        print("\n".join(set_samples(table, args.set)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
