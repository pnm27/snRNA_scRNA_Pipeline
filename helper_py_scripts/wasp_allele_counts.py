#!/usr/bin/env python3
"""Reference and alternative allele counts of a donor's reads over its
heterozygous SNPs, before and after WASP's remap test (wasp_donor.sh finish).

The SNPs and alleles are WASP's own: its snptable.py ($WASP_HOME/mapping)
reads them from the genotype h5 files, kept to those polymorphic in the
donor's genotype sample(s) (find_intersecting_snps.py's --samples), and each
read's bases over them are compared with the two alleles as
find_intersecting_snps.py does ('read SNP ref matches', ...). So the counts
are of the reads WASP tested, not of all reads:

    before  to.remap.bam: every read over a SNP that WASP remapped
    after   remap.keep.bam: those that mapped back to the same place,
            which are all the reads over SNPs in the bam after WASP (the
            reads over no SNP are merged in unchanged)

--skip_duplicates leaves out the reads marked as duplicates (0x400) from
both, as when the bam after WASP is without them (ATAC), so the change is
WASP's alone. Run with the WASP module's python (pysam, PyTables), as its
scripts are (numpy names restored as by wasp_run.py).

Usage::

    python wasp_allele_counts.py --snp_index snp_index.h5 --snp_tab snp_tab.h5 \\
        --haps haps.h5 --samples samples.txt \\
        --bam before=to.remap.bam --bam after=remap.keep.bam -o D1.allele_counts.tsv

The output is one row per bam: label, ref, alt, other (bases matching
neither allele), ref_fraction (ref / (ref + alt)) and snps (the SNPs with at
least one read).
"""

import argparse
import csv
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COLUMNS = ['label', 'ref', 'alt', 'other', 'ref_fraction', 'snps']


def read_samples(path):
    """The genotype sample IDs, one per line (find_intersecting_snps.py's --samples)."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def counted(read, skip_duplicates=False) -> bool:
    """A read whose alleles count: mapped, primary, and not a duplicate
    when those are left out."""
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary
                or (skip_duplicates and read.is_duplicate))


def count_read(read, table, counts, seen):
    """Adds a read's bases over the table's SNPs to counts [ref, alt, other]
    (as find_intersecting_snps.py's count_ref_alt_matches) and the SNPs to
    'seen' (their positions)."""
    snp_idx, read_pos, _, _ = table.get_overlapping_snps(read)
    seq = read.query_sequence
    for i, pos in zip(snp_idx, read_pos):
        base = seq[pos - 1]
        if base == table.snp_allele1[i].decode('utf-8'):
            counts[0] += 1
        elif base == table.snp_allele2[i].decode('utf-8'):
            counts[1] += 1
        else:
            counts[2] += 1
        seen.add(int(table.snp_pos[i]))


def count_bams(reads_by_label, load_table, skip_duplicates=False) -> dict:
    """{label: (ref, alt, other, snps)} of each label's reads (in the order
    of their bam, which groups them by chromosome: WASP writes to.remap.bam
    and remap.keep.bam in the order of the sorted bam). The bams are read
    side by side, a chromosome at a time, so its SNPs are loaded once
    (load_table(chrom)); a chromosome seen again is loaded again."""
    labels = list(reads_by_label)
    groups = {l: itertools.groupby(
        (r for r in reads_by_label[l] if counted(r, skip_duplicates)),
        key=lambda r: r.reference_name) for l in labels}
    heads = {l: next(groups[l], None) for l in labels}
    counts = {l: [0, 0, 0] for l in labels}
    seen = {l: {} for l in labels}
    while any(heads.values()):
        chrom = next(heads[l][0] for l in labels if heads[l])
        table = load_table(chrom)
        for l in labels:
            while heads[l] and heads[l][0] == chrom:
                snps = seen[l].setdefault(chrom, set())
                for read in heads[l][1]:
                    count_read(read, table, counts[l], snps)
                heads[l] = next(groups[l], None)
    return {l: (*counts[l], sum(len(s) for s in seen[l].values())) for l in labels}


def ref_fraction(ref, alt):
    return ref / (ref + alt) if ref + alt else None


def write(results, path):
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(COLUMNS)
        for label, (ref, alt, other, snps) in results.items():
            frac = ref_fraction(ref, alt)
            w.writerow([label, ref, alt, other, 'NA' if frac is None else f"{frac:.6g}", snps])
    os.replace(tmp, path)


def wasp_snptable(wasp_home):
    """WASP's snptable module (with the numpy names it needs restored)."""
    sys.path.insert(0, HERE)
    from wasp_run import restore_aliases
    restore_aliases()
    sys.path.insert(0, os.path.join(wasp_home, 'mapping'))
    import snptable
    return snptable


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--snp_index', required=True)
    p.add_argument('--snp_tab', required=True)
    p.add_argument('--haps', required=True)
    p.add_argument('--samples', required=True, help="file of genotype sample IDs, one per line")
    p.add_argument('--bam', action='append', required=True, metavar='LABEL=BAM',
                   help="a bam to count, with its label (before=to.remap.bam)")
    p.add_argument('--skip_duplicates', action='store_true',
                   help="leave out the reads marked as duplicates (0x400)")
    p.add_argument('--wasp_home', default=os.environ.get('WASP_HOME'))
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args(argv)
    if not args.wasp_home:
        p.error("WASP_HOME is not set: load the WASP module, or give --wasp_home")
    bams = {}
    for b in args.bam:
        label, sep, path = b.partition('=')
        if not sep or not label or not path:
            p.error(f"--bam {b}: give LABEL=BAM")
        bams[label] = path

    import pysam
    import tables
    snptable = wasp_snptable(args.wasp_home)
    samples = read_samples(args.samples)
    h5 = [tables.open_file(f, 'r') for f in (args.snp_tab, args.snp_index, args.haps)]
    table = snptable.SNPTable()

    def load_table(chrom):
        print(f"wasp_allele_counts.py: SNPs of {chrom}", file=sys.stderr)
        table.read_h5(*h5, chrom, samples)
        return table

    # WASP's bams have no index: no htslib message about it
    pysam.set_verbosity(0)
    opened = {l: pysam.AlignmentFile(p, 'rb') for l, p in bams.items()}
    try:
        results = count_bams({l: iter(f) for l, f in opened.items()}, load_table,
                             args.skip_duplicates)
    finally:
        for f in opened.values():
            f.close()
        for f in h5:
            f.close()
    write(results, args.output)
    for label, (ref, alt, other, snps) in results.items():
        frac = ref_fraction(ref, alt)
        print(f"wasp_allele_counts.py: {label}: ref {ref}, alt {alt}, other {other}, "
              f"ref fraction {'NA' if frac is None else f'{frac:.4f}'}, {snps} SNPs",
              file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
