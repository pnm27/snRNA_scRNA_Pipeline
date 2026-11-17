"""The WASP stage: samples table and paths (lib/wasp.py), stats
(helper_py_scripts/wasp_stats.py), input checks
(helper_py_scripts/wasp_check_inputs.py) and its config (schema,
schema/wasp_defaults.yaml). Needs jsonschema and pyyaml."""

import copy, csv, os, sys
import pytest
import yaml

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "helper_py_scripts"))
from lib import wasp
from lib.config_check import config_problems
from lib.config_defaults import fill_defaults, merge
import wasp_stats
import wasp_check_inputs

SCHEMA = yaml.safe_load(open(os.path.join(REPO, "schema", "config.schema.yaml")))


def table(tmp_path, text, name="samples.tsv"):
    wasp._TABLES.clear()
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# Samples table ---------------------------------------------------------------

def test_samples_table(tmp_path):
    path = table(tmp_path, "pool\tdonor\tsample\tbam\n"
                           "P1\tD1\tG1\tx/D1.bam\n"
                           "P1\tD2\tG2; G3\tx/D2.bam\n"
                           "P1\tD3\tG4,G5 ;G6\tx/D3.bam\n")
    t = wasp.read_samples_table(path, need_bam=True, sets=['main'])
    assert t[('P1', 'D1')] == {'samples': ['G1'], 'bam': 'x/D1.bam',
                               'bai': 'x/D1.bam.bai', 'genotype_set': 'main'}
    assert t[('P1', 'D2')]['samples'] == ['G2', 'G3']
    # ',' too, as WASP's own --samples string
    assert t[('P1', 'D3')]['samples'] == ['G4', 'G5', 'G6']
    assert wasp.set_samples(t, 'main') == ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']


def test_samples_table_csv_quoted_ids(tmp_path):
    path = table(tmp_path, 'pool,donor,sample\nP1,D1,"G1,G2"\nP1,D2,G3;G4\n', "samples.csv")
    t = wasp.read_samples_table(path, sets=['main'])
    assert t[('P1', 'D1')]['samples'] == ['G1', 'G2']
    assert t[('P1', 'D2')]['samples'] == ['G3', 'G4']


def test_samples_table_csv_without_bam(tmp_path):
    path = table(tmp_path, "pool,donor,sample\nP1,D1,G1\n", "samples.csv")
    t = wasp.read_samples_table(path, sets=['main'])
    assert t[('P1', 'D1')]['bam'] is None


@pytest.mark.parametrize("text, kw, error", [
    ("pool\tsample\nP1\tG1\n", {}, "missing: donor"),
    ("pool\tdonor\tsample\nP1\tD1\tG1\n", {'need_bam': True}, "missing: bam"),
    ("pool\tdonor\tsample\nP1\t\tG1\n", {}, "empty pool or donor"),
    ("pool\tdonor\tsample\nP/1\tD1\tG1\n", {}, "can't contain '/'"),
    ("pool\tdonor\tsample\nP1\tD1\tG1 G2\n", {}, "can't contain whitespace \\(G1 G2\\)"),
    ("pool\tdonor\tsample\nP1\tD1\tG1\nP1\tD1\tG2\n", {}, "listed twice: P1 D1"),
    ("pool\tdonor\tsample\nP1\tD1\tG1\n", {'sets': ['wgs', 'array']}, "no genotype_set"),
    ("pool\tdonor\tsample\tgenotype_set\nP1\tD1\tG1\tchip\n", {'sets': ['wgs']},
     "'chip' isn't one of"),
])
def test_samples_table_errors(tmp_path, text, kw, error):
    path = table(tmp_path, text)
    with pytest.raises(ValueError, match=error):
        wasp.read_samples_table(path, **kw)


def test_sample_names(tmp_path):
    names = table(tmp_path, "pool\tdonor\tsample\n"
                            "P1\t\tGENO_P1\n"        # every donor of P1
                            "P2\tD1\tGENO_A;GENO_B\n"
                            "P2\t\tGENO_P2\n", "names.tsv")
    path = table(tmp_path, "pool\tdonor\tsample\n"
                           "P1\tD1\t\n"
                           "P1\tD2\tG2\n"            # the sample column wins
                           "P2\tD1\t\n"              # (pool, donor) before pool
                           "P2\tD2\t\n"
                           "P3\tD9\t\n")             # else the donor
    t = wasp.read_samples_table(path, sets=['main'], names=names)
    assert {k: r['samples'] for k, r in t.items()} == {
        ('P1', 'D1'): ['GENO_P1'], ('P1', 'D2'): ['G2'],
        ('P2', 'D1'): ['GENO_A', 'GENO_B'], ('P2', 'D2'): ['GENO_P2'],
        ('P3', 'D9'): ['D9']}


def test_samples_table_without_sample_column(tmp_path):
    names = table(tmp_path, "pool,sample\nP1,GENO_P1\n", "names.csv")
    path = table(tmp_path, "pool\tdonor\nP1\tD1\nP2\tD2\n")
    t = wasp.read_samples_table(path, sets=['main'], names=names)
    assert t[('P1', 'D1')]['samples'] == ['GENO_P1']
    assert t[('P2', 'D2')]['samples'] == ['D2']


@pytest.mark.parametrize("text, error", [
    ("pool\tdonor\nP1\tD1\n", "missing: sample"),
    ("pool\tsample\nP1\t\n", "empty pool or sample"),
    ("pool\tsample\nP1\tG 1\n", "can't contain whitespace"),
    ("pool\tdonor\tsample\nP1\tD1\tG1\nP1\tD1\tG2\nP2\t\tG3\nP2\t\tG4\n",
     "listed twice: P1 D1, P2"),
])
def test_sample_names_errors(tmp_path, text, error):
    names = table(tmp_path, text, "names.tsv")
    with pytest.raises(ValueError, match=error):
        wasp.read_sample_names(names)


def test_set_samples_cli_with_names(tmp_path, capsys):
    names = table(tmp_path, "pool\tsample\nP1\tGENO_P1\n", "names.tsv")
    path = table(tmp_path, "pool\tdonor\nP1\tD1\nP2\tD2\n")
    wasp_check_inputs.main(['set-samples', '--table', path, '--set', 'main',
                            '--sets', 'main', '--names', names])
    assert capsys.readouterr().out.split() == ['D2', 'GENO_P1']


def test_genotype_set_column(tmp_path):
    path = table(tmp_path, "pool\tdonor\tsample\tgenotype_set\n"
                           "P1\tD1\tG1\twgs\nP1\tD2\tA2\tarray\n")
    t = wasp.read_samples_table(path, sets=['wgs', 'array'])
    assert wasp.set_samples(t, 'array') == ['A2']


def test_select_donors_and_missing_bams(tmp_path):
    (tmp_path / "D1.bam").touch()
    (tmp_path / "D1.bam.bai").touch()
    path = table(tmp_path, f"pool\tdonor\tsample\tbam\nP1\tD1\tG1\t{tmp_path}/D1.bam\n"
                           f"P1\tD2\tG2\t{tmp_path}/D2.bam\n")
    t = wasp.read_samples_table(path, need_bam=True, sets=['main'])
    assert wasp.select_donors(t, 'P1', ['D1', 'D3']) == (['D1'], ['D3'])
    assert wasp.missing_bams(t) == [f"P1 D2: no bam {tmp_path}/D2.bam",
                                    f"P1 D2: no bai {tmp_path}/D2.bam.bai"]


# Paths -------------------------------------------------------------------------

def config(**wasp_keys):
    base = {'workflow': {'start': 'align', 'aligner': 'STARsolo'},
            'modules': {'STARsolo': {'genome_pick': {'genome': 'star_idx', 'gtf': 'g.gtf'}},
                        'cellranger_arc': {'ref': {'genome': 'refs/arc/'}},
                        'wasp': {'wasp_dir': 'out', 'work_dir': 'tmp/',
                                 'genotypes': {'wgs': {'h5_dir': 'h5'},
                                               'chip': {'vcf': 'c.bcf', 'chrom_sizes': 'c.txt'}}}}}
    base['modules']['wasp'].update(wasp_keys)
    return base


def test_paths():
    c = config()
    assert wasp.donor_paths(c)['bam'] == "out/{pool}/{donor}.after_wasp.bam"
    assert wasp.work_dir(c) == "tmp/{pool}/{donor}/"
    assert wasp.genotype_files(c, 'wgs')['haps'] == "h5/haps.h5"
    assert wasp.genotype_files(c, 'chip')['snp_tab'] == "out/genotypes/chip/snp_tab.h5"
    assert wasp.vcf_sets(c) == ['chip']
    assert wasp.vcf_files(c, 'chip') == ['c.bcf']


def test_remap_reference():
    c = config()
    assert wasp.remap_reference(c) == ('star_idx', 'g.gtf')
    c['workflow']['aligner'] = 'cellranger_arc'
    assert wasp.remap_reference(c) == ('refs/arc/star', None)
    c = config(remap={'genome': 'mine', 'gtf': None})
    assert wasp.remap_reference(c) == ('mine', None)
    c = config()
    c['workflow']['start'] = 'wasp'
    with pytest.raises(ValueError, match="needs 'modules: wasp: remap: genome'"):
        wasp.remap_reference(c)


# Stats -------------------------------------------------------------------------

INTERSECT_LOG = """\
DISCARD reads:
  unmapped: 0
  mate unmapped: 0
  improper pair: 0
  different chromosome: 0
  indel: 12
  secondary alignment: 3
  supplementary alignment: 0
  excess overlapping snps: 5
  excess allelic combinations: 1
  read pairs with discordant shared SNPs: 0
  missing pairs (e.g. mismatched read names): 0
KEEP reads:
  single-end: 900
  pairs: 0
REMAP reads:
  single-end: 80
  pairs: 0
VARIANT allele counts:
  read SNP ref matches: 70
  read SNP alt matches: 30
  read SNP mismatches: 2
"""
FILTER_LOG = """\
keep reads: 72
bad reads: 6
not present reads: 0
discard reads: 2
CIGAR mismatch: 0
CIGAR missing or multiple: 1
"""


def flagstat(mapped):
    return (f"{mapped} + 0 in total (QC-passed reads + QC-failed reads)\n"
            f"{mapped} + 0 primary\n0 + 0 secondary\n"
            f"{mapped} + 0 mapped (100.00% : N/A)\n")


def logs(tmp_path):
    files = {'intersect_log': INTERSECT_LOG, 'filter_log': FILTER_LOG,
             'flagstat_before': flagstat(1000), 'flagstat_after': flagstat(972)}
    for k, v in files.items():
        (tmp_path / k).write_text(v)
    return {k: str(tmp_path / k) for k in files}


def test_parse(tmp_path):
    row = wasp_stats.parse(**logs(tmp_path))
    assert row['intersect_discard_indel'] == 12
    assert row['intersect_discard_mate_unmapped'] == 0
    assert row['intersect_keep'] == 900 and row['intersect_remap'] == 80
    assert row['filter_bad'] == 6 and row['filter_cigar_missing_or_multiple'] == 1
    assert row['summary_input'] == 1000 and row['summary_final'] == 972
    assert row['summary_failed_remap'] == 8
    assert row['summary_removed'] == 28 and row['summary_removed_pct'] == pytest.approx(2.8)
    assert row['summary_ref_allele_fraction'] == pytest.approx(0.7)


def test_parse_paired(tmp_path):
    row = wasp_stats.parse(**logs(tmp_path), paired_end=True)
    assert row['intersect_keep'] == 0  # the 'pairs:' lines
    assert row['summary_input'] == 500


def test_parse_cli_and_summarize(tmp_path):
    files = logs(tmp_path)
    for donor in ['D1', 'D2']:
        args = ['parse', '--pool', 'P1', '--donor', donor, '--sample', 'G1;G2',
                '-o', str(tmp_path / f"{donor}.tsv")]
        for k, v in files.items():
            args += [f'--{k}', v]
        assert wasp_stats.main(args) == 0
    out = tmp_path / "run"
    rows = wasp_stats.summarize(
        [str(tmp_path / "D1.tsv"), str(tmp_path / "D2.tsv")], labels=['A4', 'A13'],
        stats=str(out / "stats.tsv"), summary=str(out / "summary.tsv"))
    assert [r['dataset'] for r in rows] == ['A4', 'A13']
    with open(out / "stats.tsv") as f:
        got = list(csv.DictReader(f, delimiter='\t'))
    assert got[0]['sample'] == 'G1;G2' and got[0]['summary_final'] == '972'
    assert list(got[0])[:5] == ['dataset', 'pool', 'donor', 'sample', 'paired_end']
    with open(out / "summary.tsv") as f:
        summary = [r for r in csv.DictReader(f, delimiter='\t')
                   if r['metric'] == 'summary_final']
    assert {r['dataset'] for r in summary} == {'A4', 'A13'}
    assert summary[0]['median'] == '972' and summary[0]['sd'] == 'NA'


def test_plots(tmp_path):
    pytest.importorskip('matplotlib')
    files = logs(tmp_path)
    row = {'pool': 'P1', 'donor': 'D1', **wasp_stats.parse(**files)}
    wasp_stats.write_rows([row], list(row), str(tmp_path / "s.tsv"))
    wasp_stats.summarize([str(tmp_path / "s.tsv")], plots=str(tmp_path / "p.pdf"))
    assert (tmp_path / "p.pdf").stat().st_size > 0


# Input checks ------------------------------------------------------------------

def test_missing_samples_and_contigs(tmp_path):
    path = table(tmp_path, "pool\tdonor\tsample\nP1\tD1\tG1;G9\nP1\tD2\tG2\n")
    t = wasp.read_samples_table(path, sets=['main'])
    assert wasp_check_inputs.missing_samples(t, 'main', {'G1', 'G2'}) == ["P1 D1: G9"]
    assert wasp_check_inputs.common_contigs(['chr1', 'chr2'], {'1', '2'}) == set()
    assert wasp_check_inputs.common_contigs(['1', 'MT'], {'1', '2'}) == {'1'}


def test_check_contigs_fails(tmp_path, monkeypatch):
    (tmp_path / "contigs.txt").write_text("chr1\nchr2\n")
    monkeypatch.setattr(wasp_check_inputs, 'h5_chromosomes', lambda _: {'1', '2'})
    with pytest.raises(SystemExit, match="no chromosome in common"):
        wasp_check_inputs.check_contigs(str(tmp_path / "contigs.txt"), "snp_tab.h5", "D1.bam")


# Config ------------------------------------------------------------------------

MINIMAL = {
    'workflow': {'start': 'wasp', 'until': 'wasp'},
    'references': {'tools': {'wasp': 'wasp/0.3.4'}},
    'modules': {'wasp': {'samples': 'samples.tsv',
                         'genotypes': {'main': {'h5_dir': 'h5/'}},
                         'remap': {'genome': 'star_idx/'}}},
}


def test_minimal_wasp_config_is_valid():
    c = fill_defaults(copy.deepcopy(MINIMAL))
    assert config_problems(c, SCHEMA) == ([], [])
    assert c['workflow']['gt_check'] is False
    assert c['experiment']['demultiplexing']['method'] is None


@pytest.mark.parametrize("over, error", [
    ({'references': {'tools': {'wasp': None}}}, "WASP env module"),
    ({'modules': {'wasp': {'remap': {'genome': None}}}}, "remap: genome"),
    ({'workflow': {'until': 'split_bams'}}, "start: wasp WASP-filters"),
    ({'modules': {'wasp': {'genotypes': {'chip': {'vcf': 'x.bcf'}}}}}, "chrom_sizes"),
])
def test_wasp_config_errors(over, error):
    c = fill_defaults(merge(MINIMAL, over))
    errors, _ = config_problems(c, SCHEMA)
    assert any(error in e for e in errors), errors


def test_until_wasp_needs_gt_check_off():
    with open(os.path.join(REPO, "new_config.example.yaml")) as f:
        c = merge(yaml.safe_load(f), {
            'workflow': {'until': 'wasp', 'gt_check': True},
            'experiment': {'demultiplexing': {'method': 'vireo'}},
            'references': {'tools': {'wasp': 'wasp/0.3.4'}},
            'modules': {'wasp': {'samples': 's.tsv', 'genotypes': {'main': {'h5_dir': 'h5/'}}}}})
    errors, _ = config_problems(c, SCHEMA)
    assert any("gt_check: no" in e for e in errors), errors
    c['workflow']['gt_check'] = False
    assert not [e for e in config_problems(c, SCHEMA)[0] if 'wasp' in e.lower()]
