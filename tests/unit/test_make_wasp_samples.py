"""make_wasp_samples.py: a WASP samples table from per-donor bams (run
with pytest)."""

import os, subprocess, sys
import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
SCRIPT = os.path.join(REPO, "helper_py_scripts", "make_wasp_samples.py")
sys.path.insert(0, os.path.join(REPO, "helper_py_scripts"))
sys.path.insert(0, REPO)
from make_wasp_samples import find_bams
from lib.wasp import read_samples_table

SPLIT = "split/final"
CDNA = SPLIT + "/Sample_{pool}-cDNA/{donor}.bam"
ATAC = SPLIT + "/ATAC/Sample_{pool}-cDNA/{donor}.bam"


def touch_bams(root, paths, index=True):
    for p in paths:
        full = root / p
        full.parent.mkdir(parents=True, exist_ok=True)
        full.touch()
        if index:
            (root / f"{p}.bai").touch()


def run(root, *args):
    return subprocess.run([sys.executable, SCRIPT, *args], cwd=root,
                          capture_output=True, text=True)


@pytest.fixture
def bams(tmp_path):
    touch_bams(tmp_path, [CDNA.format(pool=p, donor=d) for p, d in
                          (("P1", "D1"), ("P1", "D2"), ("P2", "D3"))])
    touch_bams(tmp_path, [ATAC.format(pool=p, donor=d) for p, d in
                          (("P1", "D1"), ("P1", "D2"), ("P2", "D4"))])
    return tmp_path


def test_find_bams(bams, monkeypatch):
    monkeypatch.chdir(bams)
    found = find_bams(CDNA)
    assert sorted(found) == [("P1", "D1"), ("P1", "D2"), ("P2", "D3")]
    # the cDNA pattern doesn't take the ATAC bams ('/' isn't in a wildcard)
    assert all("/ATAC/" not in b for b in found.values())


def test_pattern_needs_pool_and_donor(bams):
    out = run(bams, "--cDNA", SPLIT + "/{donor}.bam")
    assert out.returncode != 0 and "needs {pool}" in out.stderr


def test_table(bams):
    out = run(bams, "--cDNA", CDNA, "--ATAC", ATAC, "-o", "wasp_samples.tsv")
    assert out.returncode == 0, out.stderr
    lines = (bams / "wasp_samples.tsv").read_text().splitlines()
    assert lines[0] == "pool\tdonor\tmodality\tbam"
    assert lines[1:3] == [f"P1\tD1\tcDNA\t{CDNA.format(pool='P1', donor='D1')}",
                          f"P1\tD1\tATAC\t{ATAC.format(pool='P1', donor='D1')}"]
    assert "cDNA bam but no ATAC bam: P2 D3" in out.stderr
    assert "ATAC bam but no cDNA bam: P2 D4" in out.stderr
    # the pipeline reads it (start: wasp)
    table = read_samples_table(str(bams / "wasp_samples.tsv"), need_bam=True)
    assert set(table) == {("P1", "D1", "cDNA"), ("P1", "D2", "cDNA"), ("P2", "D3", "cDNA"),
                          ("P1", "D1", "ATAC"), ("P1", "D2", "ATAC"), ("P2", "D4", "ATAC")}


def test_no_index_skipped(bams):
    touch_bams(bams, [CDNA.format(pool="P2", donor="D5")], index=False)
    out = run(bams, "--cDNA", CDNA, "--ATAC", ATAC, "-o", "t.tsv")
    assert "no index" in out.stderr and "D5" not in (bams / "t.tsv").read_text()
    # skipped bams aren't in the one-modality warnings
    assert "P2 D5" not in out.stderr


def test_absolute(bams):
    run(bams, "--cDNA", CDNA, "-o", "t.tsv", "--absolute")
    rows = (bams / "t.tsv").read_text().splitlines()[1:]
    assert len(rows) == 3 and all(os.path.isabs(r.split("\t")[3]) for r in rows)


def test_no_match(bams):
    out = run(bams, "--ATAC", "nowhere/{pool}/{donor}.bam")
    assert out.returncode != 0 and "No ATAC bams match" in out.stderr
