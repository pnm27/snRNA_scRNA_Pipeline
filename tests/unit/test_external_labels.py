"""Splitting bams by donor labels from outside the pipeline: the
transformation steps, normalize_external_labels.py, the external mode of
create_inp_splitBam_consolidated.py and check_external_barcodes.py."""

import os, subprocess, sys
import pandas as pd
import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
SCRIPTS = os.path.join(REPO, "helper_py_scripts")
sys.path.insert(0, SCRIPTS)
from transform_steps import apply_steps
from normalize_external_labels import normalize, read_labels
from check_external_barcodes import check

LAST_FIELD = [{'type': 'pythonic', 'operation': 'split', 'by': '_'},
              {'type': 'pythonic', 'operation': 'select_index', 'index': -1}]


# Transformation steps -------------------------------------------------------
@pytest.mark.parametrize("value, steps, out", [
    ('P12_L1_AAAC-1', LAST_FIELD, 'AAAC-1'),
    ('P12_L1_AAAC-1', [{'type': 'regex', 'pattern': '^.*_', 'replacement': ''}], 'AAAC-1'),
    ('AAAC', [{'type': 'regex', 'pattern': '$', 'replacement': '-1'}], 'AAAC-1'),
    ('AAAC-1', [{'type': 'pythonic', 'operation': 'replace', 'old': '-1', 'new': ''}], 'AAAC'),
    ('a-b-c', [{'type': 'pythonic', 'operation': 'split', 'by': '-'},
               {'type': 'pythonic', 'operation': 'select_range', 'end': -1},
               {'type': 'pythonic', 'operation': 'join', 'by': '-'}], 'a-b'),
    # Allowed by the schema, but returned None before
    (' pool ', [{'type': 'pythonic', 'operation': 'strip'},
                {'type': 'pythonic', 'operation': 'upper'}], 'POOL'),
    ('Pool', [{'type': 'regex', 'pattern': 'POOL', 'replacement': 'P',
               'flags': 'IGNORECASE'}], 'P'),
])
def test_steps(value, steps, out):
    assert apply_steps(value, steps) == out


# normalize_external_labels.py ----------------------------------------------
def labels():
    return pd.DataFrame({
        'barcodes': ['P1_L1_AAAC-1', 'P1_L1_CCCG-1', 'P1_L2_GGGT-1', 'P2_L1_AAAC-1', 'P2_L1_TTTA-1'],
        'Subj_ID': ['D1', 'D2', 'Doublet', 'D3', None],
        'pool': ['Pool1', 'Pool1', 'Pool1', 'Pool2', 'Pool2'],
    })


def test_normalize():
    out = normalize(labels(), LAST_FIELD,
                    [{'type': 'pythonic', 'operation': 'replace', 'old': 'Pool', 'new': 'P'}])
    assert out.values.tolist() == [
        ['P1', 'D1', 'AAAC-1'], ['P1', 'D2', 'CCCG-1'], ['P2', 'D3', 'AAAC-1']]


def test_normalize_same_barcode_in_a_pool():
    df = labels()
    df.loc[1, 'barcodes'] = 'P1_L2_AAAC-1'  # same cell barcode, other lane
    with pytest.raises(ValueError, match="share a barcode within a pool"):
        normalize(df, LAST_FIELD)


def test_read_table(tmp_path):
    f = tmp_path / "labels.csv.gz"
    df = labels().rename(columns={'barcodes': 'bc', 'Subj_ID': 'individualID', 'pool': 'origin'})
    df.to_csv(f, index=False)
    out = read_labels(str(f), 'individualID', 'origin', 'bc')
    assert out.columns.tolist() == ['pool', 'Subj_ID', 'barcodes'] and len(out) == 5
    with pytest.raises(ValueError, match="no column.*'nope'"):
        read_labels(str(f), 'individualID', 'nope', 'bc')
    with pytest.raises(ValueError, match="tables need 'barcode_column'"):
        read_labels(str(f), 'individualID', 'origin', None)


def test_read_h5ad(tmp_path):
    ad = pytest.importorskip("anndata")
    import numpy as np
    df = labels()
    obs = pd.DataFrame({'individualID': df['Subj_ID'].fillna('unassigned').values,
                        'origin': df['pool'].values}, index=df['barcodes'].values)
    f = tmp_path / "pooled.h5ad"
    ad.AnnData(np.zeros((5, 2), dtype=np.float32), obs=obs).write_h5ad(f)
    out = normalize(read_labels(str(f), 'individualID', 'origin', None), LAST_FIELD)
    assert out['barcodes'].tolist() == ['AAAC-1', 'CCCG-1', 'AAAC-1']


def run(script, *args, **kw):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *map(str, args)],
                          capture_output=True, text=True, **kw)


def test_normalize_cli(tmp_path):
    f, out = tmp_path / "labels.tsv", tmp_path / "out" / "external_labels.tsv"
    labels().rename(columns={'barcodes': 'barcode'}).to_csv(f, sep='\t', index=False)
    r = run("normalize_external_labels.py", f, out, '--donor_column', 'Subj_ID',
            '--pool_column', 'pool', '--barcode_column', 'barcode',
            '--barcode_steps', '[{"type": "regex", "pattern": "^.*_", "replacement": ""}]')
    assert r.returncode == 0, r.stderr
    assert pd.read_csv(out, sep='\t')['barcodes'].tolist() == ['AAAC-1', 'CCCG-1', 'AAAC-1']


# create_inp_splitBam_consolidated.py --external_pool ------------------------
def test_external_pool(tmp_path):
    inp, conv, out = tmp_path / "external_labels.tsv", tmp_path / "conv.csv", tmp_path / "P1_bc_hash.txt"
    normalize(labels(), LAST_FIELD).assign(
        pool=lambda d: d['pool'].str.replace('Pool', 'P')).to_csv(inp, sep='\t', index=False)
    pd.DataFrame({'from': ['D1'], 'to': ['R001']}).to_csv(conv, index=False)
    # the short option names lib/params.py passes
    r = run("create_inp_splitBam_consolidated.py", inp, out, '--external_pool', 'P1',
            '--conv_file', conv, '--conv_file_from_col', 'from', '--conv_file_to_col', 'to')
    assert r.returncode == 0, r.stderr
    assert pd.read_csv(out, sep='\t').values.tolist() == [['R001', 'AAAC-1'], ['D2', 'CCCG-1']]
    r = run("create_inp_splitBam_consolidated.py", inp, out, '--external_pool', 'P9')
    assert r.returncode != 0 and "No cells of pool 'P9'" in r.stderr


# check_external_barcodes.py ------------------------------------------------
def sam(*cbs):
    return ['@HD\tVN:1.6\n'] + [
        f"r{i}\t0\tchr1\t1\t255\t10M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\tCB:Z:{cb}\tUB:Z:U\n"
        for i, cb in enumerate(cbs)] + ["r_no_cb\t4\t*\t0\t0\t*\t*\t0\t0\tA\tF\tNH:i:0\n"]


def test_check_barcodes():
    ok, report = check({'AAAC-1', 'CCCG-1'}, sam('AAAC-1', 'CCCG-1', 'TTTT-1'), 0.5)
    assert ok and "2 of 3 sampled reads" in report
    # a '-1' on one side only
    ok, report = check({'AAAC', 'CCCG'}, sam('AAAC-1', 'CCCG-1'), 0.05)
    assert not ok and "AAAC-1" in report
    ok, report = check({'AAAC'}, [], 0.05)
    assert not ok and "No reads with a CB tag" in report


def test_check_barcodes_cli(tmp_path):
    h, out = tmp_path / "hash.txt", tmp_path / "check.txt"
    h.write_text("Subj_ID\tbarcodes\nD1\tAAAC\n")
    r = run("check_external_barcodes.py", '--hash_file', h, '--out', out,
            input=''.join(sam('AAAC-1')))
    assert r.returncode != 0 and "barcode_transformations" in r.stderr and not out.exists()
    r = run("check_external_barcodes.py", '--hash_file', h, '--out', out,
            input=''.join(sam('AAAC')))
    assert r.returncode == 0 and out.exists()
