"""helper_py_scripts/wasp_allele_counts.py: the ref and alt alleles of the
reads over SNPs, per bam (run with pytest; reads and WASP's SNP table are
stand-ins)."""

import os, sys
from types import SimpleNamespace as NS
import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(REPO, "helper_py_scripts"))
import wasp_allele_counts as wac


class Table:
    """As WASP's SNPTable: a read's SNPs (indices, 1-based read positions)."""
    def __init__(self, chrom):
        self.chrom = chrom
        self.snp_allele1 = [b'A', b'C']
        self.snp_allele2 = [b'G', b'T']
        self.snp_pos = [100, 200]

    def get_overlapping_snps(self, read):
        return read.snps[0], read.snps[1], [], []


def read(chrom, seq, snps, **flags):
    f = dict(is_unmapped=False, is_secondary=False, is_supplementary=False, is_duplicate=False)
    f.update(flags)
    return NS(reference_name=chrom, query_sequence=seq, snps=snps, **f)


def test_counts_per_bam_one_load_per_chromosome():
    loaded = []

    def load(chrom):
        loaded.append(chrom)
        return Table(chrom)

    before = [read('1', 'AC', ([0, 1], [1, 2])),            # ref, ref
              read('1', 'GT', ([0, 1], [1, 2])),            # alt, alt
              read('2', 'N', ([0], [1])),                   # other
              read('2', 'A', ([0], [1]), is_secondary=True),
              read('2', 'A', ([0], [1]), is_duplicate=True)]   # counted: duplicates kept
    after = [read('1', 'AC', ([0, 1], [1, 2])), read('2', 'A', ([0], [1]))]
    got = wac.count_bams({'before': iter(before), 'after': iter(after)}, load)
    # ref, alt, other, SNPs with reads
    assert got == {'before': (3, 2, 1, 3), 'after': (3, 0, 0, 3)}
    assert loaded == ['1', '2']


def test_skip_duplicates_and_chromosome_only_in_second():
    loaded = []

    def load(chrom):
        loaded.append(chrom)
        return Table(chrom)

    before = [read('1', 'A', ([0], [1])), read('1', 'G', ([0], [1]), is_duplicate=True)]
    after = [read('3', 'G', ([0], [1]))]
    got = wac.count_bams({'before': iter(before), 'after': iter(after)}, load,
                         skip_duplicates=True)
    assert got == {'before': (1, 0, 0, 1), 'after': (0, 1, 0, 1)}
    assert loaded == ['1', '3']


def test_write(tmp_path):
    out = tmp_path / "c.tsv"
    wac.write({'before': (3, 1, 0, 2), 'after': (0, 0, 0, 0)}, str(out))
    assert out.read_text() == ("label\tref\talt\tother\tref_fraction\tsnps\n"
                               "before\t3\t1\t0\t0.75\t2\nafter\t0\t0\t0\tNA\t0\n")


def test_bad_bam_argument(tmp_path):
    with pytest.raises(SystemExit):
        wac.main(['--snp_index', 'i', '--snp_tab', 't', '--haps', 'h', '--samples', 's',
                  '--bam', 'to.remap.bam', '--wasp_home', str(tmp_path), '-o', 'x'])
