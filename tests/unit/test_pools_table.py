"""lib/io.py: the pools table of 'split_by: external: pools' (the bam of
each pool, for bams outside the aligner's bams_dir)."""

import os, sys
from types import SimpleNamespace
import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, REPO)
from lib import io  # noqa: E402


def config(pools=None, input='external'):
    return {
        'workflow': {'aligner': 'STARsolo'},
        'modules': {
            'STARsolo': {'bams_dir': 'bams/', 'fold_struct': '{pool}/{pool}',
                         'bam': '_Aligned.sortedByCoord.out.bam',
                         'bai': '_Aligned.sortedByCoord.out.bam.bai'},
            'split_bams': {'temp_dir': 'tmp/', 'split_by': {
                'input': input, 'external': {'pools': pools}}},
        },
    }


@pytest.fixture(autouse=True)
def fresh_cache():
    io._POOLS_TABLES.clear()


def write(tmp_path, text, name='bams.tsv'):
    f = tmp_path / name
    f.write_text(text)
    return str(f)


def test_table_bams(tmp_path):
    path = write(tmp_path, "pool\tbam\nsamp1\t/r/dayA/samp1/a.bam\n"
                           "samp3\t/r/dayB/samp3/b.bam\n")
    c = config(path)
    wc = SimpleNamespace(pool='samp3')
    assert io.get_bam(wc, c) == {'bam': '/r/dayB/samp3/b.bam',
                                 'bai': '/r/dayB/samp3/b.bam.bai'}
    # Derived bams go to temp_dir, not next to the (maybe read-only) bam
    assert io.get_bam(wc, c, suffix='_filt.bam') == {
        'bam': 'tmp/{pool}/{pool}_filt.bam', 'bai': 'tmp/{pool}/{pool}_filt.bam.bai'}
    assert io.derived_bam_paths(c, 'cDNA', '_chr1.bam')['bam'] == 'tmp/{pool}/{pool}_chr1.bam'


def test_bai_column_and_csv(tmp_path):
    path = write(tmp_path, "pool,bam,bai\nsamp1,a.bam,idx/a.bai\nsamp2,b.bam,\n",
                 name='bams.csv')
    t = io.read_pools_table(path)
    assert t['samp1']['bai'] == 'idx/a.bai' and t['samp2']['bai'] == 'b.bam.bai'


def test_without_table_unchanged():
    for c in (config(None), config('bams.tsv', input='raw')):
        assert io.external_pools_file(c) is None
        assert io.get_bam(SimpleNamespace(pool='P1'), c, suffix='_filt.bam')['bam'] \
            == 'bams/{pool}/{pool}_filt.bam'


@pytest.mark.parametrize("text, msg", [
    ("pool\tpath\nsamp1\ta.bam\n", "missing: bam"),
    ("pool\tbam\n", "one row per pool"),
    ("pool\tbam\nsamp1\ta.bam\nsamp1\tb.bam\n", "listed twice: samp1"),
    ("pool\tbam\ndayA/samp1\ta.bam\n", "can't contain '/'"),
    ("pool\tbam\nsamp1\t\n", "empty pool or bam"),
])
def test_bad_tables(tmp_path, text, msg):
    with pytest.raises(ValueError, match=msg):
        io.read_pools_table(write(tmp_path, text))


def test_missing_pool_bams(tmp_path):
    (tmp_path / 'a.bam').write_text('')
    (tmp_path / 'a.bam.bai').write_text('')
    (tmp_path / 'b.bam').write_text('')
    path = write(tmp_path, f"pool\tbam\nsamp1\t{tmp_path}/a.bam\nsamp2\t{tmp_path}/b.bam\n")
    missing = io.missing_pool_bams(config(path), ['samp1', 'samp2', 'samp9'])
    assert missing == [f"samp2: no bai {tmp_path}/b.bam.bai",
                       f"samp9: not in the pools table {path}"]
    with pytest.raises(ValueError, match="isn't in the pools table"):
        io.pool_bam(config(path), 'samp9')
