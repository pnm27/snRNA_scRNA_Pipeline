"""rules/wasp/rna.smk, rule wasp_remap_rna's shell: a gzipped gtf is unzipped
into the work dir for STAR and removed after, also when the job fails (run
with pytest). STAR and helper_sh_scripts/wasp_donor.sh are stand-ins."""

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
"""
# prepare makes the work dir; finish removes it, as the real one does
DONOR = """#!/usr/bin/env bash
cmd=$1; shift
while [ $# -gt 0 ]; do [ "$1" = --work ] && work=$2; shift; done
[ "$cmd" = prepare ] && mkdir -p "$work"
[ "$cmd" = finish ] && rm -rf "$work"
exit 0
"""


def shell_block():
    with open(os.path.join(REPO, "rules", "wasp", "rna.smk")) as f:
        return re.search(r'shell:\s*r"""(.*?)"""', f.read(), re.S).group(1)


def run(tmp_path, gtf, fail=""):
    repo = tmp_path / "repo"
    (repo / "helper_sh_scripts").mkdir(parents=True, exist_ok=True)
    (repo / "helper_sh_scripts" / "wasp_donor.sh").write_text(DONOR)
    bin_ = tmp_path / "bin"
    bin_.mkdir(exist_ok=True)
    (bin_ / "STAR").write_text(STAR)
    (bin_ / "STAR").chmod(stat.S_IRWXU)
    work = tmp_path / "work"
    o = lambda n: str(tmp_path / n)
    cmd = snakemake_utils.format(
        shell_block(),
        params=NS(work=str(work), samples=o("s.txt"), paired="", genome=o("genome"),
                  gtf=gtf, star_extra="", fastqs="${prefix}.remap.fq.gz", threads=1,
                  donor="D1", sample_ids="G1"),
        input=NS(bam=o("D1.bam"), snp_index="i", snp_tab="t", haps="h"),
        log=NS(intersect=o("i.log"), remap=o("remap.log"), filter=o("f.log")),
        output=NS(bam=o("out.bam"), flagstat_before="fb", flagstat_after="fa", stats="st"),
        wildcards=NS(pool="P1", final_donor="D1"))
    r = subprocess.run(["bash", "-c", cmd], cwd=repo, capture_output=True, text=True,
                       env=dict(os.environ, PATH=f"{bin_}:{os.environ['PATH']}",
                                RECORD=o("star_args"), FAIL=fail))
    return r, work, (tmp_path / "star_args").read_text().split("\n", 1)


def test_gzipped_gtf_unzipped_then_removed(tmp_path):
    gz = tmp_path / "genes.gtf.gz"
    with gzip.open(gz, "wt") as f:
        f.write(GTF)
    r, work, (given, text) = run(tmp_path, str(gz))
    assert r.returncode == 0, r.stderr
    assert given == f"{work}/remap.gtf" and text == GTF
    assert not os.path.exists(given)


def test_gzipped_gtf_removed_when_star_fails(tmp_path):
    gz = tmp_path / "genes.gtf.gz"
    with gzip.open(gz, "wt") as f:
        f.write(GTF)
    r, work, (given, _) = run(tmp_path, str(gz), fail="star")
    assert r.returncode != 0
    assert work.is_dir() and not (work / "remap.gtf").exists()   # work dir kept, gtf not


def test_plain_gtf_as_given(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)
    r, work, (given, text) = run(tmp_path, str(gtf))
    assert r.returncode == 0, r.stderr
    assert given == str(gtf) and text == GTF and gtf.exists()


def test_no_gtf(tmp_path):
    r, work, (given, _) = run(tmp_path, "")
    assert r.returncode == 0, r.stderr
    assert given == ""
