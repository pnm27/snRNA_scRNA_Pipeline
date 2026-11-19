#!/usr/bin/env python3
"""Write a WASP samples table (modules: wasp: samples) from per-donor bams.

Each modality gets a path pattern with {pool} and {donor}; every bam that
matches it becomes a row (pool, donor, modality, bam), e.g. the outputs of
split_bams:

  python helper_py_scripts/make_wasp_samples.py \\
      --cDNA 'projdir/STARsolo_split_filt_bams/final/Sample_{pool}-cDNA/{donor}.bam' \\
      --ATAC 'projdir/STARsolo_split_filt_bams/final/ATAC/Sample_{pool}-cDNA/{donor}.bam' \\
      -o wasp_samples.tsv

A bam with no index (<bam>.bai) is skipped with a warning. Donors found in
only one modality are reported (they get WASP of that modality only).
"""
import argparse
import glob
import os
import re
import sys


def find_bams(pattern: str) -> dict:
    """{(pool, donor): bam} of the files matching 'pattern'."""
    for key in ('{pool}', '{donor}'):
        if key not in pattern:
            sys.exit(f"Pattern {pattern} needs {key}")
    parts = re.split(r'(\{pool\}|\{donor\})', pattern)
    regex = ''.join(
        {'{pool}': r'(?P<pool>[^/]+)', '{donor}': r'(?P<donor>[^/]+)'}.get(p, re.escape(p))
        for p in parts)
    gpattern = pattern.replace('{pool}', '*').replace('{donor}', '*')
    found = {}
    for path in sorted(glob.glob(gpattern)):
        m = re.fullmatch(regex, path)
        if m:
            found[(m['pool'], m['donor'])] = path
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cDNA', help="pattern of the cDNA bams, with {pool} and {donor}")
    ap.add_argument('--ATAC', help="pattern of the ATAC bams, with {pool} and {donor}")
    ap.add_argument('-o', '--out', default='wasp_samples.tsv')
    ap.add_argument('--absolute', action='store_true', help="write absolute bam paths")
    args = ap.parse_args()
    patterns = {m: p for m, p in (('cDNA', args.cDNA), ('ATAC', args.ATAC)) if p}
    if not patterns:
        sys.exit("Give --cDNA and/or --ATAC")

    rows, per_modality = [], {}
    for modality, pattern in patterns.items():
        bams = find_bams(pattern)
        if not bams:
            sys.exit(f"No {modality} bams match {pattern}")
        per_modality[modality] = set()
        for (pool, donor), bam in sorted(bams.items()):
            if not os.path.exists(f"{bam}.bai"):
                print(f"WARNING: no index {bam}.bai: skipped", file=sys.stderr)
                continue
            per_modality[modality].add((pool, donor))
            rows.append((pool, donor, modality,
                         os.path.abspath(bam) if args.absolute else bam))

    if len(per_modality) == 2:
        for m, other in (('cDNA', 'ATAC'), ('ATAC', 'cDNA')):
            only = sorted(per_modality[m] - per_modality[other])
            if only:
                print(f"WARNING: {m} bam but no {other} bam: "
                      + ", ".join(f"{p} {d}" for p, d in only), file=sys.stderr)

    rows.sort(key=lambda r: (r[0], r[1], r[2] != 'cDNA'))
    with open(args.out, 'w') as f:
        f.write("pool\tdonor\tmodality\tbam\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{args.out}: {len(rows)} bams, "
          f"{len({(r[0], r[1]) for r in rows})} donors, "
          f"{len({r[0] for r in rows})} pools", file=sys.stderr)


if __name__ == '__main__':
    main()
