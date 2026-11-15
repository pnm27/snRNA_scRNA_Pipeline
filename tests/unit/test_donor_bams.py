"""helper_py_scripts/donor_bams.py: per-donor bams named by the swap
correction file, linked or merged (samtools is faked: merge writes the
list of its inputs, index touches the .bai, view -H prints the file)."""

import os
import stat
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'helper_py_scripts'))

import donor_bams  # noqa: E402

SQ = "@SQ\tSN:chr1\tLN:100\n"

FAKE_SAMTOOLS = """#!/usr/bin/env bash
case "$1" in
    view) cat "$3" ;;
    index) shift; while [[ "$1" == -* ]]; do shift 2; done; touch "$1.bai" ;;
    merge) shift; out=""; ins=()
        while [[ $# -gt 0 ]]; do
            case "$1" in -o) out=$2; shift 2 ;; -@) shift 2 ;; -*) shift ;;
                *) ins+=("$(basename "$1")"); shift ;; esac
        done
        printf "%s\\n" "${ins[@]}" > "$out" ;;
esac
"""


@pytest.fixture
def project(tmp_path, monkeypatch):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    fake = bin_dir / 'samtools'
    fake.write_text(FAKE_SAMTOOLS)
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv('PATH', f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.chdir(tmp_path)

    for pool, donor in [('P1', 'v1'), ('P1', 'v2'), ('P2', 'v1'), ('P2', 'v3'),
                        ('P3', 'v4')]:
        d = tmp_path / 'split' / f"Sample_{pool}-cDNA"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{donor}.bam").write_text(SQ)
        (d / f"{donor}.bam.bai").write_text('')

    pd.DataFrame({
        'pool':   ['P1', 'P1', 'P2', 'P2', 'P3'],
        'vcf_id': ['v1', 'v2', 'v1', 'v3', 'v4'],
        'indiv':  ['A', 'B', 'A', '', 'C'],
        'status': ['True', 'yes', 'True', 'True', 'False'],
        'tissue': ['DLPFC', 'DLPFC', 'AMY', 'AMY', 'AMY'],
    }).to_csv(tmp_path / 'swap.csv', index=False)
    return tmp_path


def args(*extra):
    return ['--swap', 'swap.csv', '--pool_column', 'pool',
            '--donor_column', 'vcf_id', '--name_column', 'indiv',
            '--status_column', 'status',
            '--bams', 'split/Sample_{pool}-cDNA/{donor}.bam',
            '--out_dir', 'out', *extra]


def manifest():
    return pd.read_csv('out/manifest.tsv', sep='\t', dtype=str, keep_default_na=False)


def test_merge_by_name_and_drop(project):
    donor_bams.main(args())
    # A: P1 and P2 merged, each with its pool's read group
    assert open('out/A.bam').read().split() == ['P1.bam', 'P2.bam']
    assert os.path.exists('out/A.bam.bai')
    assert not os.path.exists('out/A.bam.tmp_inputs')
    # B: one bam, hard linked with its index
    assert os.path.samefile('out/B.bam', 'split/Sample_P1-cDNA/v2.bam')
    assert os.path.samefile('out/B.bam.bai', 'split/Sample_P1-cDNA/v2.bam.bai')
    # empty name (P2 v3) and status False (P3 v4) dropped
    assert sorted(os.listdir('out')) == ['A.bam', 'A.bam.bai', 'B.bam', 'B.bam.bai',
                                          'manifest.tsv']
    m = manifest()
    assert dict(zip(m['output'], m['action'])) == {
        'out/A.bam': 'merge', 'out/B.bam': 'link', '': 'dropped'}
    assert (m['action'] == 'dropped').sum() == 2


def test_no_combine(project):
    donor_bams.main(args('--no_combine'))
    assert os.path.samefile('out/A/P1.bam', 'split/Sample_P1-cDNA/v1.bam')
    assert os.path.samefile('out/A/P2.bam', 'split/Sample_P2-cDNA/v1.bam')
    assert os.path.samefile('out/B/P1.bam', 'split/Sample_P1-cDNA/v2.bam')


def test_keep_apart(project):
    donor_bams.main(args('--keep_apart', 'tissue'))
    assert os.path.samefile('out/A_DLPFC.bam', 'split/Sample_P1-cDNA/v1.bam')
    assert os.path.samefile('out/A_AMY.bam', 'split/Sample_P2-cDNA/v1.bam')
    assert os.path.exists('out/B_DLPFC.bam')


def test_dry_run_writes_nothing(project, capsys):
    donor_bams.main(args('--dry_run'))
    assert not os.path.exists('out')
    assert 'out/A.bam\tmerge' in capsys.readouterr().out


def test_missing_bams_listed(project):
    os.remove('split/Sample_P1-cDNA/v2.bam')
    with pytest.raises(SystemExit, match='v2.bam'):
        donor_bams.main(args())


def test_different_references_not_merged(project):
    open('split/Sample_P2-cDNA/v1.bam', 'w').write("@SQ\tSN:1\tLN:100\n")
    with pytest.raises(SystemExit, match='different references'):
        donor_bams.main(args())


def test_existing_kept_unless_force(project):
    donor_bams.main(args())
    open('out/A.bam', 'w').write('old')
    donor_bams.main(args())
    assert open('out/A.bam').read() == 'old'
    donor_bams.main(args('--force'))
    assert open('out/A.bam').read().split() == ['P1.bam', 'P2.bam']


def test_two_donors_of_a_pool_one_name(project):
    df = pd.read_csv('swap.csv', dtype=str, keep_default_na=False)
    df.loc[1, 'indiv'] = 'A'
    df.to_csv('swap.csv', index=False)
    donor_bams.main(args())
    assert open('out/A.bam').read().split() == ['P1_v1.bam', 'P1_v2.bam', 'P2.bam']
    donor_bams.main(args('--no_combine', '--force'))
    assert sorted(os.listdir('out/A')) == ['P1_v1.bam', 'P1_v1.bam.bai',
                                           'P1_v2.bam', 'P1_v2.bam.bai',
                                           'P2.bam', 'P2.bam.bai']


def test_missing_column(project):
    with pytest.raises(SystemExit, match='no column'):
        donor_bams.main(args('--keep_apart', 'region'))
