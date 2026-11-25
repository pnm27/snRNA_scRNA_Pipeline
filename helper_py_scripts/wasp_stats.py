#!/usr/bin/env python3
"""Statistics of WASP filtering (stage 'wasp', rules/wasp/).

parse      one row for a donor bam, from the logs of
           find_intersecting_snps.py and filter_remapped_reads.py and
           samtools flagstat of the bam before and after WASP
           (python standard library only: runs in the WASP module)
summarize  joins rows (per-donor files, or the stats tables of whole
           runs, e.g. of several datasets; --collect: every donor's
           file in a modality's folder, from all batches run into it),
           and writes the table, a
           summary (mean, median, sd, min, max of each summary_* column,
           per dataset) and plots (histograms, smoothed by a kernel
           density, of the summary_* columns that tell donors apart:
           PLOT_COLUMNS; then the reference allele fraction before and
           after WASP, see below)

Counts are reads for single-end bams (sc/snRNA-seq: only R2 has the
cDNA) and pairs for paired-end ones (ATAC). With --flagstat_wasp (the bam
after WASP, before the duplicates were removed: ATAC), summary_wasp_kept
and summary_duplicates_removed tell WASP's part from the duplicates'.

The reference allele fraction (ref / (ref + alt) of the reads' bases over
the donor's heterozygous SNPs):

    summary_ref_allele_fraction          every read over a SNP in the bam
                                         before WASP (find_intersecting_snps.py's
                                         log)
    summary_ref_allele_fraction_tested   the reads WASP remapped, before its
                                         test; with --allele_counts
                                         (wasp_allele_counts.py)
    summary_ref_allele_fraction_final    those of them it kept: all the reads
                                         over SNPs in the bam after WASP
    summary_ref_allele_fraction_change   final - tested
    summary_ref_removed_pct,             the % of the ref and of the alt
    summary_alt_removed_pct              bases that WASP removed (a ref bias
                                         removes more ref than alt)

Usage::

    python3 wasp_stats.py parse --pool P1 --donor D1 --final_donor D1 --modality cDNA \\
        --sample GENO1 \\
        --intersect_log D1.find_intersecting_snps.log \\
        --filter_log D1.filter_remapped_reads.log \\
        --flagstat_before D1.before_wasp.flagstat \\
        --flagstat_after D1.after_wasp.flagstat -o D1.wasp_stats.tsv

    # every donor of a modality's folder, of all batches (rule wasp_stats)
    python3 wasp_stats.py summarize --collect wasp_bams/cDNA -o wasp_bams/cDNA

    # several datasets, one label each
    python3 wasp_stats.py summarize A4/wasp_stats.tsv A13/wasp_stats.tsv \\
        --labels A4 A13 -o wasp_all/
"""

import argparse
import csv
import glob
import math
import os
import re
import statistics
import sys
import time


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


def allele_counts(path) -> dict:
    """{label: {ref, alt, other, snps}} of wasp_allele_counts.py's table."""
    with open(path, newline='') as f:
        return {r['label']: {k: _int(r.get(k)) for k in ('ref', 'alt', 'other', 'snps')}
                for r in csv.DictReader(f, delimiter='\t')}


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def allele_summary(counts) -> tuple:
    """The alleles_* columns and the summary ones of the reference allele
    fraction before (the reads WASP tested) and after its test."""
    row, summary = {}, {}
    for label in ('before', 'after'):
        for k, v in counts.get(label, {}).items():
            row[f'alleles_{label}_{k}'] = v
    b, a = counts.get('before', {}), counts.get('after', {})
    tested = _div(b.get('ref'), _sum(b.get('ref'), b.get('alt')))
    final = _div(a.get('ref'), _sum(a.get('ref'), a.get('alt')))
    summary['ref_allele_fraction_tested'] = tested
    summary['ref_allele_fraction_final'] = final
    summary['ref_allele_fraction_change'] = (
        None if tested is None or final is None else final - tested)
    for k in ('ref', 'alt'):
        kept = _div(a.get(k), b.get(k))
        summary[f'{k}_removed_pct'] = None if kept is None else 100 * (1 - kept)
    return row, summary


def parse(intersect_log, filter_log, flagstat_before, flagstat_after,
          paired_end=False, flagstat_wasp=None, allele_counts_path=None) -> dict:
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
    if allele_counts_path:
        counts_row, counts_summary = allele_summary(allele_counts(allele_counts_path))
        row.update(counts_row)
        summary.update(counts_summary)
    row.update({f'summary_{k}': v for k, v in summary.items()})
    return row


def _fmt(v):
    if v is None:
        return 'NA'
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def _tmp(path):
    """A temporary name next to 'path' (same folder, same extension), to
    write to and then move over 'path': a reader (another batch's run) never
    sees a half-written file."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    stem, ext = os.path.splitext(path)
    return f"{stem}.tmp{os.getpid()}{ext}"


def write_rows(rows, columns, path):
    tmp = _tmp(path)
    with open(tmp, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r.get(c)) for c in columns])
    os.replace(tmp, path)


# A donor's stats file under a modality's folder: <pool>/stats/<final_donor>.wasp_stats.tsv
DONOR_STATS = os.path.join('*', 'stats', '*.wasp_stats.tsv')


def collect(dirs, inputs=()):
    """'inputs' and every donor's stats file under each of 'dirs' (a
    modality's folder in wasp_dir: all batches run into it), each once,
    sorted by path (pool, then final donor)."""
    found = {os.path.abspath(p): p for p in inputs}
    for d in dirs:
        for p in glob.glob(os.path.join(d, DONOR_STATS)):
            found.setdefault(os.path.abspath(p), p)
    return [found[k] for k in sorted(found)]


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


def kde(xs, grid):
    """Gaussian kernel density of xs at the grid points (bandwidth by
    Scott's rule), or None for fewer than 2 values or all the same."""
    import numpy as np
    xs = np.asarray(xs, dtype=float)
    if len(xs) < 2 or xs.std(ddof=1) == 0:
        return None
    bw = xs.std(ddof=1) * len(xs) ** (-1 / 5)
    z = (np.asarray(grid)[:, None] - xs[None, :]) / bw
    return np.exp(-0.5 * z ** 2).sum(axis=1) / (len(xs) * bw * math.sqrt(2 * math.pi))


def _values(rows, col, ds):
    return [x for r in rows if (r.get('dataset') or 'all') == ds
            for x in [_number(r.get(col))] if x is not None]


def _hist_kde(ax, rows, col, datasets, colors):
    """Histograms of a column per dataset (the same bins), each with its
    kernel density scaled to the counts."""
    import numpy as np
    values = {ds: _values(rows, col, ds) for ds in datasets}
    every = [x for xs in values.values() for x in xs]
    if not every:
        ax.text(0.5, 0.5, 'no values', ha='center', va='center', transform=ax.transAxes)
        return
    lo, hi = min(every), max(every)
    if lo == hi:
        lo, hi = lo - 0.5, hi + 0.5
    edges = np.linspace(lo, hi, 31)
    pad = (hi - lo) * 0.1
    grid = np.linspace(lo - pad, hi + pad, 300)
    for ds in datasets:
        xs = values[ds]
        if not xs:
            continue
        ax.hist(xs, bins=edges, alpha=0.3, color=colors[ds], label=f"{ds} (n={len(xs)})")
        dens = kde(xs, grid)
        if dens is not None:
            ax.plot(grid, dens * len(xs) * (edges[1] - edges[0]), color=colors[ds], lw=1.5)


def _ref_fraction_page(fig, rows, datasets, colors):
    """The reference allele fraction of the reads WASP tested, before and
    after its test (per donor bam), and the share of ref and alt bases it
    removed; False if no row has them."""
    import numpy as np
    pairs = {ds: [(b, a) for r in rows if (r.get('dataset') or 'all') == ds
                  for b, a in [(_number(r.get('summary_ref_allele_fraction_tested')),
                                _number(r.get('summary_ref_allele_fraction_final')))]
                  if b is not None and a is not None]
             for ds in datasets}
    if not any(pairs.values()):
        return False
    every = [x for ps in pairs.values() for p in ps for x in p] + [0.5]
    lo, hi = min(every), max(every)
    pad = max((hi - lo) * 0.15, 0.01)
    grid = np.linspace(lo - pad, hi + pad, 400)
    axes = fig.subplots(2, 2)
    fig.suptitle("Reference allele fraction over the donor's heterozygous SNPs, "
                 "of the reads WASP tested", fontsize=11, fontweight='bold')

    ax = axes[0, 0]
    for ds in datasets:
        if not pairs[ds]:
            continue
        for i, (name, ls) in enumerate((('before WASP', '--'), ('after WASP', '-'))):
            xs = [p[i] for p in pairs[ds]]
            dens = kde(xs, grid)
            if dens is not None:
                ax.plot(grid, dens, color=colors[ds], ls=ls, lw=1.6, label=f"{ds}: {name}")
            else:
                ax.axvline(xs[0], color=colors[ds], ls=ls, label=f"{ds}: {name}")
    ax.axvline(0.5, color='grey', lw=0.8, ls=':')
    ax.set_title('per donor bam (kernel density)', fontsize=10)
    ax.set_xlabel('ref / (ref + alt)')
    ax.set_ylabel('density')
    ax.legend(fontsize=7)

    ax = axes[0, 1]
    for ds in datasets:
        if pairs[ds]:
            b, a = zip(*pairs[ds])
            ax.scatter(b, a, s=10, alpha=0.6, color=colors[ds], label=ds)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color='grey', lw=0.8, ls=':')
    ax.axhline(0.5, color='grey', lw=0.5, ls=':')
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_title('each donor bam: before vs after', fontsize=10)
    ax.set_xlabel('before WASP (reads tested)')
    ax.set_ylabel('after WASP (reads kept)')

    ax = axes[1, 0]
    _hist_kde(ax, rows, 'summary_ref_allele_fraction_change', datasets, colors)
    ax.axvline(0, color='grey', lw=0.8, ls=':')
    ax.set_title('change: after - before', fontsize=10)
    ax.set_xlabel('change in ref / (ref + alt)')
    ax.set_ylabel('donor bams')

    ax = axes[1, 1]
    removed = []
    for ds in datasets:
        pts = [(x, y) for r in rows if (r.get('dataset') or 'all') == ds
               for x, y in [(_number(r.get('summary_alt_removed_pct')),
                             _number(r.get('summary_ref_removed_pct')))]
               if x is not None and y is not None]
        if pts:
            removed += [v for p in pts for v in p]
            x, y = zip(*pts)
            ax.scatter(x, y, s=10, alpha=0.6, color=colors[ds], label=ds)
    if removed:
        top = max(removed) * 1.05 or 1
        ax.plot([0, top], [0, top], color='grey', lw=0.8, ls=':')
        ax.set_xlim(0, top)
        ax.set_ylim(0, top)
    ax.set_title('bases over SNPs removed by WASP', fontsize=10)
    ax.set_xlabel('alt bases removed (%)')
    ax.set_ylabel('ref bases removed (%)')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return True


# The first page of the plots: the summary_* columns that tell donors apart
# (the others are counts of these, or their complements: final_pct =
# 100 - removed_pct; the ref fraction before and after WASP has its own page)
PLOT_COLUMNS = {
    'summary_input': 'reads before WASP (ATAC: pairs)',
    'summary_sent_for_remapping_pct': '% over a SNP: remapped and tested',
    'summary_failed_remap_pct': '% failed the remap test',
    'summary_removed_pct': '% removed in all (final_pct = 100 - this)',
    'summary_duplicates_removed_pct': '% duplicates removed (ATAC)',
    'summary_ref_allele_fraction': 'ref allele fraction, all reads over SNPs, before WASP',
}


def plot(rows, path):
    """A pdf: histograms of the PLOT_COLUMNS that have values, one colour
    per dataset, each smoothed by its kernel density; then, if the rows
    have them, the reference allele fraction before and after WASP."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    cols = [c for c in PLOT_COLUMNS if any(_number(r.get(c)) is not None for r in rows)]
    datasets = sorted({r.get('dataset') or 'all' for r in rows})
    palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = {ds: palette[i % len(palette)] for i, ds in enumerate(datasets)}
    tmp = _tmp(path)
    with PdfPages(tmp) as pdf:
        ncol = 2
        nrow = max(1, math.ceil(len(cols) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(12, 3.2 * nrow), squeeze=False)
        for ax, col in zip(axes.flat, cols):
            _hist_kde(ax, rows, col, datasets, colors)
            ax.set_title(col.replace('summary_', ''), fontsize=10, fontweight='bold')
            ax.set_xlabel(PLOT_COLUMNS[col], fontsize=8)
            ax.set_ylabel('donor bams')
        if 'summary_ref_allele_fraction' in cols:
            ax = axes.flat[cols.index('summary_ref_allele_fraction')]
            ax.axvline(0.5, color='grey', lw=0.8, ls=':')
        for ax in list(axes.flat)[len(cols):]:
            ax.axis('off')
        handles, labels = axes.flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc='upper center', ncol=min(len(labels), 4))
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        pdf.savefig(fig)
        plt.close(fig)

        fig = plt.figure(figsize=(12, 9))
        if _ref_fraction_page(fig, rows, datasets, colors):
            pdf.savefig(fig)
        plt.close(fig)
    os.replace(tmp, path)


def summarize(inputs, labels=None, stats=None, summary=None, plots=None):
    if not inputs:
        raise ValueError("no stats files to summarize")
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


def _time(t=None):
    return time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t))


def write_log(path, inputs, labels, rows, out):
    """Appends to 'path' the stats files joined (each with its time and
    label), the number of rows and the files written, so a run's outputs
    can be traced to what they were made from."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'a') as f:
        f.write(f"[{_time()}] wasp_stats.py summarize: {len(inputs)} stats file(s), "
                f"{len(rows)} row(s)\n")
        f.write("stats files used (modified, path"
                + (", label" if labels else "") + "):\n")
        for i, p in enumerate(inputs):
            f.write(f"  {_time(os.path.getmtime(p))}\t{os.path.abspath(p)}"
                    + (f"\t{labels[i]}" if labels else "") + "\n")
        f.write("written:\n")
        for k, p in out.items():
            if p:
                f.write(f"  {k}\t{os.path.abspath(p)}\n")


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
    a.add_argument('--allele_counts',
                   help="wasp_allele_counts.py's table (before: the reads WASP tested, after: kept)")
    a.add_argument('-o', '--output', required=True)

    s = sub.add_parser('summarize', help="join stats files; summary and plots")
    s.add_argument('inputs', nargs='*', help="stats files (per donor, or of whole runs)")
    s.add_argument('--collect', nargs='+', metavar='DIR',
                   help="also every donor's stats file under DIR "
                        "(<pool>/stats/*.wasp_stats.tsv: a modality's folder in wasp_dir, "
                        "so the donors of every batch run into it)")
    s.add_argument('--labels', nargs='+', help="a dataset label for each input")
    s.add_argument('-o', '--out_dir', help="writes wasp_stats.tsv, wasp_summary.tsv, wasp_plots.pdf")
    s.add_argument('--stats')
    s.add_argument('--summary')
    s.add_argument('--plots')
    s.add_argument('--log', help="appends what was joined: the time, each stats file "
                                 "(with its time) and the files written")
    args = p.parse_args(argv)

    if args.cmd == 'parse':
        row = {'pool': args.pool, 'donor': args.donor,
               'final_donor': args.final_donor or args.donor, 'modality': args.modality,
               'sample': args.sample, 'paired_end': args.paired_end}
        row.update(parse(args.intersect_log, args.filter_log, args.flagstat_before,
                         args.flagstat_after, args.paired_end, args.flagstat_wasp,
                         args.allele_counts))
        write_rows([row], list(row), args.output)
        return 0

    out = {k: getattr(args, k) for k in ('stats', 'summary', 'plots')}
    if args.out_dir:
        names = {'stats': 'wasp_stats.tsv', 'summary': 'wasp_summary.tsv', 'plots': 'wasp_plots.pdf'}
        out = {k: v or os.path.join(args.out_dir, names[k]) for k, v in out.items()}
    if not any(out.values()):
        p.error("give -o/--out_dir, or --stats, --summary and/or --plots")
    inputs = args.inputs
    if args.collect:
        if args.labels:
            p.error("--labels label the inputs one by one: not with --collect")
        inputs = collect(args.collect, inputs)
    if not inputs:
        p.error("no stats files: give them, or --collect a folder that has them")
    print(f"wasp_stats.py summarize: {len(inputs)} stats file(s)", file=sys.stderr)
    rows = summarize(inputs, args.labels, **out)
    if args.log:
        write_log(args.log, inputs, args.labels, rows, out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
