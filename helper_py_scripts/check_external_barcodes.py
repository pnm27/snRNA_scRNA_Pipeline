#!/usr/bin/env python3
"""Check that the barcodes of a donor hash file are the bam's CB tags.

For bams split by donor labels from outside the pipeline ('split_by:
external'): the barcodes there were transformed back from unique names
(e.g. 'P12_L1_AAACCTGAGC-1'), so a wrong transformation, or a '-1'
present on only one side (cellranger vs STARsolo), would give empty
per-donor bams after hours of work. This stops the run early instead.

Reads SAM lines (e.g. the first reads of 'samtools view') from stdin and
counts those whose CB tag is a barcode of the hash file. Fails when fewer
than --min_fraction of the reads with a CB tag match.

Usage::

    (samtools view pool.bam || true) | head -n 200000 | \\
        python3 check_external_barcodes.py --hash_file P1_bc_hash.txt \\
            --out P1_bc_check.txt
"""

from __future__ import annotations
import argparse, sys
from typing import Iterable, TextIO


def read_hash_barcodes(path: str) -> set[str]:
    """Barcodes (2nd column) of a donor hash file (with a header)."""
    with open(path) as f:
        next(f, None)
        return {line.rstrip('\n').split('\t')[1] for line in f if line.strip()}


def cb_tags(sam_lines: Iterable[str]) -> Iterable[str]:
    """The CB tag of each read that has one."""
    for line in sam_lines:
        if line.startswith('@'):
            continue
        for field in line.rstrip('\n').split('\t')[11:]:
            if field.startswith('CB:Z:'):
                yield field[5:]
                break


def check(barcodes: set[str], sam_lines: Iterable[str],
          min_fraction: float) -> tuple[bool, str]:
    """(passed, report) of the reads' CB tags against 'barcodes'."""
    n = matched = 0
    unmatched: list[str] = []
    for cb in cb_tags(sam_lines):
        n += 1
        if cb in barcodes:
            matched += 1
        elif len(unmatched) < 3 and cb not in unmatched and cb != '-':
            unmatched.append(cb)
    frac = matched / n if n else 0.0
    report = (
        f"{matched} of {n} sampled reads with a CB tag ({frac:.1%}) are "
        f"from cells of the hash file (needed: {min_fraction:.0%}).\n"
        f"  bam CB tags, e.g.:        {', '.join(unmatched) or '-'}\n"
        f"  hash file barcodes, e.g.: {', '.join(sorted(barcodes)[:3])}"
    )
    if n == 0:
        return False, ("No reads with a CB tag in the sample (is the bam "
                       "empty, or without CB tags?)\n" + report)
    return frac >= min_fraction, report


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="Check that a hash file's barcodes are the bam's CB tags "
        "(SAM lines on stdin).")
    parser.add_argument('--hash_file', required=True)
    parser.add_argument('--out', required=True, help="Report, written on success")
    parser.add_argument('--min_fraction', type=float, default=0.05,
                        help="Least fraction of sampled reads (with a CB tag) "
                        "from the hash file's cells (default: 0.05). Reads of "
                        "empty droplets also carry CB tags, so it's low; a "
                        "barcode format mismatch gives ~0")
    return parser


def main(stdin: TextIO = sys.stdin):
    args = get_argument_parser().parse_args()
    passed, report = check(read_hash_barcodes(args.hash_file), stdin,
                           args.min_fraction)
    print(report)
    if not passed:
        sys.exit(
            f"The barcodes of {args.hash_file} don't look like the bam's CB "
            "tags. Check 'modules: split_bams: split_by: external: "
            "barcode_transformations' (e.g. a '-1' suffix on one side only).")
    with open(args.out, 'w') as f:
        f.write(report + '\n')


if __name__ == '__main__':
    main()
