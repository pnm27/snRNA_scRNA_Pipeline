#!/usr/bin/env python3
"""Checks of the WASP stage's inputs (rules/wasp/common.smk), before the long jobs.

samples
    each genotype sample ID of the samples table is in its set's
    haps.h5 (else find_intersecting_snps.py fails hours in)
contigs
    a bam and a set's snp_tab.h5 have chromosomes in common, as WASP
    matches them: the same name, or with the 'chr' prefix added or
    stripped (chr1 vs 1); with none in common, WASP would find no SNP
    and keep every read
set-samples
    the genotype sample IDs of a set (not run by the pipeline: e.g. to
    subset VCFs by hand with bcftools view -S)

The h5 files are read with PyTables (or h5py), as in the WASP module's
python. snp2h5 writes a table 'samples_<chrom>' (column 'name') per
chromosome in haps.h5 and a node per chromosome in snp_tab.h5.

Usage::

    python wasp_check_inputs.py samples --table samples.tsv \\
        --set main=/geno/haps.h5,/geno/snp_tab.h5 [--names names.tsv] -o checked.txt
    python wasp_check_inputs.py contigs --contigs contigs.txt \\
        --snp_tab /geno/snp_tab.h5 --bam donor.bam
    python wasp_check_inputs.py set-samples --table samples.tsv --set main --sets main \\
        [--names names.tsv]

--names is the sample names table ('modules: wasp: sample_names'), for the
rows of the samples table with no sample.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lib.wasp import read_samples_table, set_samples, row_name  # noqa: E402


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
    """'<pool> <donor>[ <modality>]: <sample>' of the set's samples not in 'known'."""
    return [f"{row_name(key)}: {s}"
            for key, row in sorted(table.items())
            if row['genotype_set'] == gset
            for s in row['samples'] if s not in known]


def h5_chrom_of(contig: str, h5_chroms) -> str:
    """The h5 chromosome WASP reads for a bam contig (snptable.read_h5): the
    same name, else with 'chr' stripped (chr1 -> 1) or added (1 -> chr1);
    None if neither is there."""
    if contig in h5_chroms:
        return contig
    other = contig[3:] if contig.startswith('chr') else f"chr{contig}"
    return other if other in h5_chroms else None


def common_contigs(bam_contigs, h5_chroms) -> set:
    """The bam contigs WASP finds SNPs of."""
    return {c for c in bam_contigs if h5_chrom_of(c, h5_chroms)}


def check_samples(table_path, sets: dict, need_bam: bool, out=None, names=None):
    table = read_samples_table(table_path, need_bam=need_bam, sets=list(sets), names=names)
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
            f"WASP samples table {table_path}: the 'sample' column (or the "
            "sample names table, or else the donor name) must give the "
            "genotype IDs of each set (find_intersecting_snps.py's "
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
            "WASP matches names as they are or with 'chr' added or stripped.")
    renamed = sorted(c for c in common if c not in chroms)
    if renamed:
        print(f"{bam}: {len(renamed)} chromosome(s) matched to {snp_tab} by the 'chr' "
              f"prefix, as WASP does (e.g. {renamed[0]} -> {h5_chrom_of(renamed[0], chroms)})",
              file=sys.stderr)
    return common


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('samples')
    s.add_argument('--table', required=True)
    s.add_argument('--set', action='append', required=True,
                   help="<name>=<haps.h5>,<snp_tab.h5>, once per genotype set")
    s.add_argument('--need_bam', action='store_true')
    s.add_argument('--names', help="the sample names table")
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
    g.add_argument('--names', help="the sample names table")
    args = p.parse_args(argv)

    if args.cmd == 'samples':
        sets = {}
        for item in args.set:
            name, files = item.split('=', 1)
            haps, snp_tab = files.split(',', 1)
            sets[name] = (haps, snp_tab)
        check_samples(args.table, sets, args.need_bam, args.output, args.names)
    elif args.cmd == 'contigs':
        check_contigs(args.contigs, args.snp_tab, args.bam)
    else:
        table = read_samples_table(args.table, need_bam=args.need_bam, sets=args.sets,
                                   names=args.names)
        print("\n".join(set_samples(table, args.set)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
