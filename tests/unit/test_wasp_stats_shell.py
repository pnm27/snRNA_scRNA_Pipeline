"""rules/wasp/common.smk, rule wasp_stats's shell (run with pytest; needs
matplotlib): each run's log is named by its start time and lists the stats
files used and the files written; a failure names the log."""

import os, re, subprocess, sys
import pytest

pytest.importorskip("matplotlib")

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ROW = ("pool\tdonor\tsummary_input\tsummary_removed_pct\n"
       "{pool}\t{donor}\t{n}\t{pct}\n")


def shell_block():
    with open(os.path.join(REPO, "rules", "wasp", "common.smk")) as f:
        text = f.read()
    rule = text[text.index("rule wasp_stats:"):]
    return re.search(r'shell:\s*r"""(.*?)"""', rule, re.S).group(1)


def run(tmp_path, folders):
    bin_ = tmp_path / "bin"
    bin_.mkdir(exist_ok=True)
    (bin_ / "python3").symlink_to(sys.executable)
    # as Snakemake fills the params in (the shell has no other braces)
    params = {'folders': " ".join(str(f) + "/" for f in folders),
              'log_dir': f"{tmp_path}/wasp/logs/", 'configfiles': "/c/batch1.yaml"}
    cmd = shell_block()
    for k, v in params.items():
        cmd = cmd.replace("{params.%s}" % k, v)
    assert "{" not in cmd.replace("${", "")
    return subprocess.run(["bash", "-c", cmd], cwd=REPO, capture_output=True, text=True,
                          env=dict(os.environ, PATH=f"{bin_}:{os.environ['PATH']}"))


def donor(folder, pool, d, n):
    p = folder / pool / "stats" / f"{d}.wasp_stats.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(ROW.format(pool=pool, donor=d, n=n, pct=n / 100))
    return p


def test_log_lists_inputs_and_outputs(tmp_path):
    cdna = tmp_path / "wasp" / "cDNA"
    files = [donor(cdna, "P1", "D1", 1000), donor(cdna, "P2", "D3", 1500)]
    r = run(tmp_path, [cdna])
    assert r.returncode == 0, r.stderr
    logs = list((tmp_path / "wasp" / "logs").glob("wasp_stats.*.log"))
    assert len(logs) == 1 and re.fullmatch(r"wasp_stats\.\d{8}_\d{6}\.log", logs[0].name)
    assert f"wasp_stats: log {logs[0]}" in r.stderr
    log = logs[0].read_text()
    assert "wasp_stats started" in log and "wasp_stats done" in log
    assert "config: /c/batch1.yaml" in log and "commit " in log
    assert "2 stats file(s), 2 row(s)" in log
    for f in files:
        assert f"\t{f}\n" in log
    for name in ("wasp_stats.tsv", "wasp_summary.tsv", "wasp_plots.pdf"):
        assert (cdna / name).exists() and f"\t{cdna / name}\n" in log


def test_failure_names_the_log(tmp_path):
    cdna = tmp_path / "wasp" / "cDNA"
    cdna.mkdir(parents=True)   # no stats files: summarize fails
    r = run(tmp_path, [cdna])
    assert r.returncode != 0
    assert "wasp_stats: failed" in r.stderr and "wasp_stats." in r.stderr
