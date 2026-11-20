"""helper_py_scripts/wasp_run.py: WASP's scripts run with numpy's removed
aliases (np.int) restored (run with pytest)."""

import os, subprocess, sys
import pytest

np = pytest.importorskip("numpy")
REPO = os.path.join(os.path.dirname(__file__), "..", "..")
WRAPPER = os.path.join(REPO, "helper_py_scripts", "wasp_run.py")

# As WASP's mapping/: a script importing a module of its folder that uses np.int
SNPTABLE = "import numpy as np\ndef idx(xs):\n    return np.array(xs, dtype=np.int)\n"
SCRIPT = """import sys, snptable
print(sys.argv[1:], snptable.idx([1, 2]).dtype.kind, __name__)
"""
sys.path.insert(0, os.path.join(REPO, "helper_py_scripts"))
import wasp_run


@pytest.fixture
def mapping(tmp_path):
    d = tmp_path / "WASP" / "mapping"
    d.mkdir(parents=True)
    (d / "snptable.py").write_text(SNPTABLE)
    (d / "find_intersecting_snps.py").write_text(SCRIPT)
    return d


def run(*args, cwd):
    return subprocess.run([sys.executable, *args], cwd=cwd, capture_output=True, text=True)


def test_runs_script_with_np_int(mapping, tmp_path):
    out = run(WRAPPER, str(mapping / "find_intersecting_snps.py"), "--is_sorted", "a.bam",
              cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "['--is_sorted', 'a.bam'] i __main__"


@pytest.mark.skipif(hasattr(np, "int"), reason="numpy < 1.24 still has np.int")
def test_fails_without_wrapper(mapping, tmp_path):
    # what the WASP module's python does with numpy >= 1.24
    out = run(str(mapping / "find_intersecting_snps.py"), cwd=tmp_path)
    assert out.returncode != 0 and "has no attribute 'int'" in out.stderr


def test_script_failure_passes_through(mapping, tmp_path):
    (mapping / "bad.py").write_text("import sys\nsys.exit(3)\n")
    assert run(WRAPPER, str(mapping / "bad.py"), cwd=tmp_path).returncode == 3


def test_no_script(tmp_path):
    out = run(WRAPPER, "nowhere.py", cwd=tmp_path)
    assert out.returncode != 0 and "no script nowhere.py" in out.stderr


def test_in1d_as_numpy_had_it():
    a, b = np.array([[1, 2], [3, 4]]), [2, 4, 9]
    assert wasp_run.in1d(a, b).tolist() == [False, True, False, True]  # flattened
    assert wasp_run.in1d([1, 5], b, invert=True).tolist() == [True, True]


def test_restore_only_missing(monkeypatch):
    monkeypatch.delattr(np, "in1d", raising=False)
    restored = wasp_run.restore_aliases()
    assert "in1d" in restored and np.in1d is wasp_run.in1d
    assert np.in1d(np.array([3, 7]), [7]).tolist() == [False, True]
    # names numpy still has are left alone
    assert "int32" not in restored and "isin" not in restored


def test_no_numpy_warnings_in_logs(mapping, tmp_path):
    # np.in1d is deprecated (numpy 2.0-2.3) or gone (2.4): either way, no
    # warning from the wrapper or the script lands in WASP's log
    (mapping / "uses_in1d.py").write_text(
        "import numpy as np\nprint(np.in1d([1, 2], [2]).tolist(), np.int(3))\n")
    out = run(WRAPPER, str(mapping / "uses_in1d.py"), cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[False, True] 3" and out.stderr == ""
