"""rules/wasp/rna.smk, rule wasp_remap_rna's shell (run with pytest; STAR is a
stand-in): STAR's bam is the output and its other files are removed; a
gzipped gtf is unzipped into the scratch dir for STAR and removed after,
also when STAR fails (the scratch dir is then kept)."""

import gzip, os, re, stat, subprocess
from types import SimpleNamespace as NS
import pytest

snakemake_utils = pytest.importorskip("snakemake.utils")

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
GTF = 'chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "G1";\n'

# Records the gtf STAR is given and what it holds; FAIL=star fails it
STAR = """#!/usr/bin/env bash
gtf=""; out=""
while [ $# -gt 0 ]; do
  case $1 in --sjdbGTFfile) gtf=$2 ;; --outFileNamePrefix) out=$2 ;; esac; shift
done
echo "$gtf" > "$RECORD"
[ -n "$gtf" ] && cat "$gtf" >> "$RECORD"
[ "$FAIL" = star ] && exit 1
echo "Uniquely mapped reads | 10" > "${out}Log.final.out"
echo bam > "${out}Aligned.sortedByCoord.out.bam"
echo log > "${out}Log.out"
"""


def shell_block():
    with open(os.path.join(REPO, "rules", "wasp", "rna.smk")) as f:
        return re.search(r'shell:\s*r"""(.*?)"""', f.read(), re.S).group(1)


def run(tmp_path, gtf, fail=""):
    bin_ = tmp_path / "bin"
    bin_.mkdir(exist_ok=True)
    (bin_ / "STAR").write_text(STAR)
    (bin_ / "STAR").chmod(stat.S_IRWXU)
    scratch = tmp_path / "work" / "remap"
    # a leftover of an earlier try, cleared first
    scratch.mkdir(parents=True, exist_ok=True)
    (scratch / "old").write_text("x")
    o = lambda n: str(tmp_path / n)
    cmd = snakemake_utils.format(
        shell_block(),
        params=NS(scratch=f"{scratch}/", genome=o("genome"), gtf=gtf, star_extra="",
                  threads=1),
        input=NS(fastqs=o("remap.fq.gz")),
        log=o("remap.log"),
        output=NS(bam=o("work/remapped.bam")),
        wildcards=NS(pool="P1", final_donor="D1"))
    r = subprocess.run(["bash", "-c", cmd], cwd=tmp_path, capture_output=True, text=True,
                       env=dict(os.environ, PATH=f"{bin_}:{os.environ['PATH']}",
                                RECORD=o("star_args"), FAIL=fail))
    return r, scratch, (tmp_path / "star_args").read_text().split("\n", 1)


def test_gzipped_gtf_unzipped_then_removed(tmp_path):
    gz = tmp_path / "genes.gtf.gz"
    with gzip.open(gz, "wt") as f:
        f.write(GTF)
    r, scratch, (given, text) = run(tmp_path, str(gz))
    assert r.returncode == 0, r.stderr
    assert given == f"{scratch}/remap.gtf" and text == GTF
    assert not os.path.exists(given)


def test_gzipped_gtf_removed_when_star_fails(tmp_path):
    gz = tmp_path / "genes.gtf.gz"
    with gzip.open(gz, "wt") as f:
        f.write(GTF)
    r, scratch, (given, _) = run(tmp_path, str(gz), fail="star")
    assert r.returncode != 0
    # the scratch dir kept (a fresh one), the gtf not
    assert scratch.is_dir() and not (scratch / "old").exists()
    assert not (scratch / "remap.gtf").exists()
    assert not (tmp_path / "work" / "remapped.bam").exists()


def test_plain_gtf_as_given(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)
    r, scratch, (given, text) = run(tmp_path, str(gtf))
    assert r.returncode == 0, r.stderr
    assert given == str(gtf) and text == GTF and gtf.exists()


def test_no_gtf(tmp_path):
    r, scratch, (given, _) = run(tmp_path, "")
    assert r.returncode == 0, r.stderr
    assert given == ""


def test_output_bam_and_scratch_removed(tmp_path):
    r, scratch, _ = run(tmp_path, "")
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "work" / "remapped.bam").read_text() == "bam\n"
    assert not scratch.exists()
    log = (tmp_path / "remap.log").read_text()
    assert "Uniquely mapped reads | 10" in log
