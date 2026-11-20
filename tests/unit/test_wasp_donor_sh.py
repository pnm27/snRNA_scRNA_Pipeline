"""helper_sh_scripts/wasp_donor.sh: each step's output and failures go to its
log (run with pytest). python and samtools are stand-ins on PATH."""

import os, stat, subprocess
import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
SCRIPT = os.path.join(REPO, "helper_sh_scripts", "wasp_donor.sh")

# FAIL names the stand-in that fails: contigs, intersect or filter
PYTHON = """#!/usr/bin/env bash
case "$*" in
  *wasp_check_inputs.py*)
    echo "D1.bam: 2 chromosome(s) matched by the 'chr' prefix" >&2
    [[ $FAIL == contigs ]] && { echo "no chromosome in common" >&2; exit 1; }; exit 0 ;;
  *find_intersecting_snps.py*)
    echo "KEEP reads:"
    [[ $FAIL == intersect ]] && { echo "Traceback: KeyError" >&2; exit 1; }; exit 0 ;;
  *filter_remapped_reads.py*)
    echo "keep reads: 10"
    [[ $FAIL == filter ]] && { echo "Traceback: IOError" >&2; exit 1; }; exit 0 ;;
esac
"""
SAMTOOLS = """#!/usr/bin/env bash
[[ $1 == view ]] && printf '@SQ\\tSN:1\\tLN:10\\n'
[[ $1 == merge && $FAIL == merge ]] && { echo "[merge] truncated file" >&2; exit 1; }
exit 0
"""


@pytest.fixture
def env(tmp_path):
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    for name, text in (("python", PYTHON), ("samtools", SAMTOOLS)):
        (bin_ / name).write_text(text)
        (bin_ / name).chmod(stat.S_IRWXU)
    return dict(os.environ, PATH=f"{bin_}:{os.environ['PATH']}", WASP_HOME=str(tmp_path))


def prepare(tmp_path, env, fail=""):
    return subprocess.run(
        ["bash", SCRIPT, "prepare", "--work", str(tmp_path / "work"), "--bam", "D1.bam",
         "--samples", "G1", "--snp_index", "i.h5", "--snp_tab", "t.h5", "--haps", "h.h5",
         "--log", str(tmp_path / "intersect.log")],
        env=dict(env, FAIL=fail), capture_output=True, text=True)


def finish(tmp_path, env, fail=""):
    (tmp_path / "work").mkdir(exist_ok=True)
    return subprocess.run(
        ["bash", SCRIPT, "finish", "--work", str(tmp_path / "work"), "--bam", "D1.bam",
         "--remapped", "r.bam", "--out", str(tmp_path / "out.bam"), "--threads", "1",
         "--filter_log", str(tmp_path / "filter.log"), "--intersect_log", "i.log",
         "--flagstat_before", str(tmp_path / "fb"), "--flagstat_after", str(tmp_path / "fa"),
         "--pool", "P1", "--donor", "D1", "--sample", "G1", "--modality", "cDNA",
         "--stats", str(tmp_path / "stats.tsv")],
        env=dict(env, FAIL=fail), capture_output=True, text=True)


def test_prepare_logs(tmp_path, env):
    out = prepare(tmp_path, env)
    assert out.returncode == 0, out.stderr
    log = (tmp_path / "intersect.log").read_text()
    # the contigs check's note and find_intersecting_snps.py's output
    assert "matched by the 'chr' prefix" in log and "KEEP reads:" in log
    assert out.stderr == ""


@pytest.mark.parametrize("fail, message", [
    ("contigs", "no chromosome in common"),
    ("intersect", "Traceback: KeyError"),
])
def test_prepare_failure_in_log(tmp_path, env, fail, message):
    out = prepare(tmp_path, env, fail)
    assert out.returncode != 0
    log = (tmp_path / "intersect.log").read_text()
    assert message in log and "wasp_donor.sh prepare: failed (exit 1)" in log
    assert "wasp_donor.sh prepare: failed" in out.stderr and "intersect.log" in out.stderr


@pytest.mark.parametrize("fail, message", [
    ("filter", "Traceback: IOError"),
    ("merge", "[merge] truncated file"),
])
def test_finish_failure_in_log(tmp_path, env, fail, message):
    out = finish(tmp_path, env, fail)
    assert out.returncode != 0
    log = (tmp_path / "filter.log").read_text()
    assert "keep reads: 10" in log and message in log
    assert "wasp_donor.sh finish: failed" in log and "filter.log" in out.stderr
    # the work dir is kept when a step fails
    assert (tmp_path / "work").is_dir()
