#!/usr/bin/env python3
"""Statistics of WASP filtering (stage 'wasp', rules/wasp/).

parse      one row for a donor bam, from the logs of
           find_intersecting_snps.py and filter_remapped_reads.py and
           samtools flagstat of the bam before and after WASP
           (python standard library only: runs in the WASP module)
summarize  joins rows (per-donor files, or the stats tables of whole
           runs, e.g. of several datasets), and writes the table, a
           summary (mean, median, sd, min, max of each summary_* column,
           per dataset) and plots (histograms of the summary_* columns)

Counts are reads for single-end bams (sc/snRNA-seq: only R2 has the
cDNA) and pairs for paired-end ones (ATAC). With --flagstat_wasp (the bam
after WASP, before the duplicates were removed: ATAC), summary_wasp_kept
and summary_duplicates_removed tell WASP's part from the duplicates'.

Usage::

    python3 wasp_stats.py parse --pool P1 --donor D1 --final_donor D1 --modality cDNA \\
        --sample GENO1 \\
        --intersect_log D1.find_intersecting_snps.log \\
        --filter_log D1.filter_remapped_reads.log \\
        --flagstat_before D1.before_wasp.flagstat \\
        --flagstat_after D1.after_wasp.flagstat -o D1.wasp_stats.tsv

    # several datasets, one label each
    python3 wasp_stats.py summarize A4/wasp_stats.tsv A13/wasp_stats.tsv \\
        --labels A4 A13 -o wasp_all/
"""

import argparse
import csv
import math
import os
import re
import statistics
import sys


# find_intersecting_snps.py's log: column -> label of its line
INTERSECT_LINES = {
    'discard_unmapped': 'unmapped',
    'discard_mate_unmapped': 'mate unmapped',
    'discard_improper_pair': 'improper pair',
    'discard_different_chromosome': 'different chromosome',
    'discard_indel': 'indel',
    'discard_secondary': 'secondary alignment',
    'discard_supplementary': 'supplementary alignment',
    'discard_excess_snps': 'excess overlapping snps',
    'discard_excess_allelic_combinations': 'excess allelic combinations',
    'discard_discordant_shared_snps': 'read pairs with discordant shared SNPs',
    'discard_missing_pairs': 'missing pairs (e.g. mismatched read names)',
    'snp_ref_matches': 'read SNP ref matches',
    'snp_alt_matches': 'read SNP alt matches',
    'snp_mismatches': 'read SNP mismatches',
}
# Blocks with a 'single-end:' and a 'pairs:' line
INTERSECT_BLOCKS = {'keep': 'KEEP reads:', 'remap': 'REMAP reads:'}
# filter_remapped_reads.py's log
FILTER_LINES = {
    'keep': 'keep reads',
    'bad': 'bad reads',
    'not_present': 'not present reads',
    'discard': 'discard reads',
    'cigar_mismatch': 'CIGAR mismatch',
    'cigar_missing_or_multiple': 'CIGAR missing or multiple',
}

ID_COLUMNS = ['dataset', 'pool', 'donor', 'final_donor', 'modality', 'sample', 'paired_end']


def single_line(lines, label):
    """The number on the first line starting with 'label' (optional ':')."""
    pat = re.compile(r'^\s*' + re.escape(label) + r'\s*:?\s*([0-9,]+)\s*$')
    for line in lines:
        m = pat.match(line)
        if m:
            return int(m.group(1).replace(',', ''))
    return None


def block_line(lines, header, kind):
    """The '<kind>:' number in the lines after 'header' (e.g. KEEP reads:)."""
    pat = re.compile(r'^\s*' + re.escape(kind) + r'\s*:\s*([0-9,]+)')
    for i, line in enumerate(lines):
        if line.strip() == header:
            for nxt in lines[i + 1:i + 7]:
                m = pat.match(nxt)
                if m:
                    return int(m.group(1).replace(',', ''))
            return None
    return None


def flagstat(lines) -> dict:
    """total and mapped (QC-passed) of samtools flagstat."""
    out = {'total': None, 'mapped': None}
    for line in lines:
        m = re.match(r'^\s*(\d+) \+ \d+ (.*)$', line)
        if not m:
            continue
        n, what = int(m.group(1)), m.group(2)
        if out['total'] is None and what.startswith('in total'):
            out['total'] = n
        elif out['mapped'] is None and re.match(r'mapped \(', what):
            out['mapped'] = n
    return out


def _read(path):
    with open(path) as f:
        return f.read().splitlines()


def _div(a, b):
    return None if a is None or not b else a / b


def _sum(*xs):
    return None if any(x is None for x in xs) else sum(xs)


def parse(intersect_log, filter_log, flagstat_before, flagstat_after,
          paired_end=False, flagstat_wasp=None) -> dict:
    """The metrics of one bam (see the module's docstring)."""
    il, fl = _read(intersect_log), _read(filter_log)
    before, after = flagstat(_read(flagstat_before)), flagstat(_read(flagstat_after))
    kind = 'pairs' if paired_end else 'single-end'
    unit = 2 if paired_end else 1

    row = {f'intersect_{k}': single_line(il, label) for k, label in INTERSECT_LINES.items()}
    for k, header in INTERSECT_BLOCKS.items():
        row[f'intersect_{k}'] = block_line(il, header, kind)
    row.update({f'filter_{k}': single_line(fl, label) for k, label in FILTER_LINES.items()})
    row.update({
        'flagstat_before_total': before['total'],
        'flagstat_before_mapped': before['mapped'],
        'flagstat_after_total': after['total'],
        'flagstat_after_mapped': after['mapped'],
    })

    n_in = _div(before['mapped'], unit)
    n_out = _div(after['mapped'], unit)
    failed = _sum(row['filter_bad'], row['filter_discard'])
    ref, alt = row['intersect_snp_ref_matches'], row['intersect_snp_alt_matches']
    summary = {
        'input': n_in,
        'no_snp_kept': row['intersect_keep'],
        'sent_for_remapping': row['intersect_remap'],
        'failed_remap': failed,
        'final': n_out,
        'removed': None if n_in is None or n_out is None else n_in - n_out,
    }
    keys = ['no_snp_kept', 'sent_for_remapping', 'failed_remap', 'final', 'removed']
    if flagstat_wasp:
        wasp = flagstat(_read(flagstat_wasp))
        row['flagstat_wasp_mapped'] = wasp['mapped']
        n_wasp = _div(wasp['mapped'], unit)
        summary['wasp_kept'] = n_wasp
        summary['duplicates_removed'] = (
            None if n_wasp is None or n_out is None else n_wasp - n_out)
        keys += ['wasp_kept', 'duplicates_removed']
    for k in keys:
        frac = _div(summary[k], n_in)
        summary[f'{k}_pct'] = None if frac is None else 100 * frac
    summary['ref_allele_fraction'] = _div(ref, _sum(ref, alt))
    row.update({f'summary_{k}': v for k, v in summary.items()})
    return row


def _fmt(v):
    if v is None:
        return 'NA'
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_rows(rows, columns, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r.get(c)) for c in columns])


def read_rows(path, label=None):
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if label is not None:
        for r in rows:
            r['dataset'] = label
    return rows


def _number(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def summary_columns(rows):
    cols = []
    for r in rows:
        cols += [c for c in r if c.startswith('summary_') and c not in cols]
    return cols


def summarize_rows(rows) -> list:
    """mean, median, sd, min and max of each summary_* column, per dataset."""
    out = []
    datasets = sorted({r.get('dataset') or 'all' for r in rows})
    for col in summary_columns(rows):
        for ds in datasets:
            xs = [x for r in rows if (r.get('dataset') or 'all') == ds
                  for x in [_number(r.get(col))] if x is not None]
            out.append({
                'dataset': ds, 'metric': col, 'n': len(xs),
                'mean': statistics.fmean(xs) if xs else None,
                'median': statistics.median(xs) if xs else None,
                'sd': statistics.stdev(xs) if len(xs) > 1 else None,
                'min': min(xs) if xs else None,
                'max': max(xs) if xs else None,
            })
    return out


def plot(rows, path):
    """Histograms of the summary_* columns, one colour per dataset."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    cols = summary_columns(rows)
    datasets = sorted({r.get('dataset') or 'all' for r in rows})
    ncol = 2
    nrow = max(1, math.ceil(len(cols) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.6 * nrow), squeeze=False)
    for ax, col in zip(axes.flat, cols):
        for ds in datasets:
            xs = [x for r in rows if (r.get('dataset') or 'all') == ds
                  for x in [_number(r.get(col))] if x is not None]
            if xs:
                ax.hist(xs, bins=30, alpha=0.4, label=f"{ds} (n={len(xs)})")
        ax.set_title(col.replace('summary_', ''), fontsize=10, fontweight='bold')
        ax.set_ylabel('donor bams')
    for ax in list(axes.flat)[len(cols):]:
        ax.axis('off')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc='upper center', ncol=min(len(labels), 4))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def summarize(inputs, labels=None, stats=None, summary=None, plots=None):
    if labels and len(labels) != len(inputs):
        raise ValueError(f"{len(labels)} labels for {len(inputs)} stats files")
    rows = []
    for i, path in enumerate(inputs):
        rows += read_rows(path, labels[i] if labels else None)
    columns = [c for c in ID_COLUMNS if any(c in r for r in rows)]
    for r in rows:
        columns += [c for c in r if c not in columns]
    if stats:
        write_rows(rows, columns, stats)
    if summary:
        write_rows(summarize_rows(rows),
                   ['dataset', 'metric', 'n', 'mean', 'median', 'sd', 'min', 'max'], summary)
    if plots:
        plot(rows, plots)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)

    a = sub.add_parser('parse', help="one row for a donor bam")
    a.add_argument('--pool', required=True)
    a.add_argument('--donor', required=True, help="the donor as split")
    a.add_argument('--final_donor', help="the donor after swap correction (default: --donor)")
    a.add_argument('--modality', default='cDNA')
    a.add_argument('--sample', required=True, help="genotype sample ID(s), ';'-separated")
    a.add_argument('--paired_end', action='store_true')
    for k in ['intersect_log', 'filter_log', 'flagstat_before', 'flagstat_after']:
        a.add_argument(f'--{k}', required=True)
    a.add_argument('--flagstat_wasp', help="flagstat after WASP, before removing duplicates")
    a.add_argument('-o', '--output', required=True)

    s = sub.add_parser('summarize', help="join stats files; summary and plots")
    s.add_argument('inputs', nargs='+', help="stats files (per donor, or of whole runs)")
    s.add_argument('--labels', nargs='+', help="a dataset label for each input")
    s.add_argument('-o', '--out_dir', help="writes wasp_stats.tsv, wasp_summary.tsv, wasp_plots.pdf")
    s.add_argument('--stats')
    s.add_argument('--summary')
    s.add_argument('--plots')
    args = p.parse_args(argv)

    if args.cmd == 'parse':
        row = {'pool': args.pool, 'donor': args.donor,
               'final_donor': args.final_donor or args.donor, 'modality': args.modality,
               'sample': args.sample, 'paired_end': args.paired_end}
        row.update(parse(args.intersect_log, args.filter_log, args.flagstat_before,
                         args.flagstat_after, args.paired_end, args.flagstat_wasp))
        write_rows([row], list(row), args.output)
        return 0

    out = {k: getattr(args, k) for k in ('stats', 'summary', 'plots')}
    if args.out_dir:
        names = {'stats': 'wasp_stats.tsv', 'summary': 'wasp_summary.tsv', 'plots': 'wasp_plots.pdf'}
        out = {k: v or os.path.join(args.out_dir, names[k]) for k, v in out.items()}
    if not any(out.values()):
        p.error("give -o/--out_dir, or --stats, --summary and/or --plots")
    summarize(args.inputs, args.labels, **out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
