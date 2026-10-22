#!/usr/bin/env python3
"""Convert a config from the old (flat) layout to the new one.

New layout:
  references  genomes, barcodes, regions and tool versions
  workflow    aligner, until (replaces 'last_step'), picard_metrics, gt_check
  experiment  inputs, hto, demultiplexing (replaces 'demux_type',
              'demux_run_type' and 'hto_demux_type'), qc
  modules     the '*_pipeline' blocks, each with its own 'fold_struct*'
  execution   rule resources (was 'resources')

Comments inside the moved blocks and the anchors (e.g. '*3v3') are kept.
Configs that only override some keys (e.g. the dry-run scenarios) are
converted too: only the keys they have are moved.

Usage:
python3 migrate_config.py new_config.yaml -o new_config.migrated.yaml
python3 migrate_config.py tests/dryrun/scenarios/*.yaml --in-place
"""

import argparse, io, sys
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.constructor import DuplicateKeyError


# Old key (top level) -> new path. Keys handled by 'derived' are not here.
MOVES = {
    # references
    'barcodes':                  'references.barcodes',
    'genome_builds':             'references.genome_builds',
    'STARsolo_genome_dir':       'references.STARsolo_genome_dir',
    'cellranger_arc_genome_dir': 'references.cellranger_arc_genome_dir',
    'cellranger_genome_dir':     'references.cellranger_genome_dir',
    'whitelist':                 'references.whitelist',
    'reg_chr_bed':               'references.reg_chr_bed',
    'chr_prefix':                'references.chr_prefix',
    'mito_prefix':               'references.mito_prefix',
    'STAR_version':              'references.tools.star',
    'STAR_snv_aware_version':    'references.tools.star_snv_aware',
    'cellranger_version':        'references.tools.cellranger',
    'cellranger-arc_version':    'references.tools.cellranger_arc',
    'featuremap_script':         'references.tools.featuremap_script',
    # workflow
    'picard_metrics':            'workflow.picard_metrics',
    'gt_check':                  'workflow.gt_check',
    # experiment
    'select_fastqs':             'experiment.inputs.select_fastqs',
    'wildcards_select':          'experiment.inputs.wildcards_select',
    'wet_lab_info':              'experiment.inputs.wet_lab_info',
    'cDNA_fastqs_dir':           'experiment.inputs.fastqs_dir.cDNA',
    'HTO_fastqs_dir':            'experiment.inputs.fastqs_dir.HTO',
    'ATAC_fastqs_dir':           'experiment.inputs.fastqs_dir.ATAC',
    'R1_suffix':                 'experiment.inputs.R1_suffix',
    'R2_suffix':                 'experiment.inputs.R2_suffix',
    'multiome_suffix':           'experiment.inputs.multiome_suffix',
    'run_tag':                   'experiment.demultiplexing.run_tag',
    'create_final_h5ad':         'experiment.demultiplexing.create_final_h5ad',
    'max_mito_percentage':       'experiment.qc.max_mito_percentage',
    'min_genes_per_cell':        'experiment.qc.min_genes_per_cell',
    'min_cells_per_gene':        'experiment.qc.min_cells_per_gene',
    # modules (a block's keys move with it)
    'STARsolo_pipeline':         'modules.STARsolo',
    'cellranger_count':          'modules.cellranger',
    'cellranger_arc_count':      'modules.cellranger_arc',
    'picard_pipeline':           'modules.picard',
    'kb_pipeline':               'modules.kite',
    'hashsolo_demux_pipeline':   'modules.hashsolo',
    'gt_demux_pipeline':         'modules.gt_demux',
    'demultiplex':               'modules.demultiplex',
    'split_bams_pipeline':       'modules.split_bams',
    'identify_swaps':            'modules.identify_swaps',
    'multiBamSummary':           'modules.deeptools.multiBamSummary',
    'plotCorrelation':           'modules.deeptools.plotCorrelation',
    # folder structures, in the module that writes them
    'fold_struct':               'modules.STARsolo.fold_struct',
    'fold_struct_kb':            'modules.kite.fold_struct',
    'fold_struct_demux':         'modules.hashsolo.fold_struct',
    'fold_struct_gt_demux':      'modules.gt_demux.fold_struct',
    'fold_struct_filt_bc':       'modules.gt_demux.fold_struct_filt_bc',
    'fold_struct_final':         'modules.demultiplex.fold_struct',
    'fold_struct_bam_split1':    'modules.split_bams.fold_struct_hash',
    'fold_struct_bam_split2':    'modules.split_bams.fold_struct_bams',
    'fold_struct_gt_demux_redo': 'modules.split_bams.fold_struct_donor',
    'fold_struct_swaps_check':   'modules.identify_swaps.fold_struct',
    'fold_struct_deeptools':     'modules.deeptools.fold_struct',
    # execution
    'resources':                 'execution.resources',
}

# Old keys no code reads: dropped
UNUSED = {'fold_struct_gt_demux2'}

SECTIONS = {
    'references': "Genomes, barcodes, regions and tool versions",
    'workflow':   "What to run: the aligner and the last stage ('until')",
    'experiment': "What was sequenced and how it's multiplexed",
    'modules':    "Parameters, output dirs and folder structures of each module",
    'execution':  "Rule resources (override rules/resources.smk)",
}

# Order of the keys made by the migration (other keys follow)
ORDER = {
    'workflow':                 ['aligner', 'until', 'picard_metrics', 'gt_check'],
    'experiment':               ['inputs', 'hto', 'demultiplexing', 'qc'],
    'experiment.demultiplexing': ['method', 'add_to_previous', 'run_tag',
                                  'genotype', 'create_final_h5ad'],
    'experiment.demultiplexing.genotype': ['enabled', 'multi_vcf'],
    'experiment.hto':           ['enabled', 'multiplexing'],
}

# The stages 'until' can name, in order
STAGES = ['align', 'cellsnp', 'demux', 'split_bams', 'identify_swaps']

# last_step values whose stage doesn't follow from their name
SPECIAL_STEPS = {
    'starsolo_kb_solo':            'align',
    'starsolo_cellsnp':            'cellsnp',
    'starsolo_demux':              'demux',
    'starsolo_gt_demux_multi_vcf': 'demux',
    'all_multi_vcf':               'demux',
    'cellranger_cellsnp':          'cellsnp',
    'cellranger_demux':            'demux',
    'multiome_alignment':          'align',
    'multiome_gt_demux':           'demux',
}


def aligner(step: str) -> str:
    """Aligner of a last_step (its prefix; 'all_multi_vcf' is STARsolo)."""
    step = step.lower()
    if step.startswith('multiome'):
        return 'cellranger_arc'
    if step.startswith('cellranger'):
        return 'cellranger'
    if step.startswith('starsolo') or step == 'all_multi_vcf':
        return 'STARsolo'
    raise ValueError(f"Unrecognised last_step: '{step}'")


def until(step: str) -> str:
    """Last stage of a last_step (as rules/produce_targets.smk reads it)."""
    step = step.lower()
    aligner(step)
    if step in SPECIAL_STEPS:
        return SPECIAL_STEPS[step]
    for stage in ['identify_swaps', 'split_bams']:
        if stage in step:
            return stage
    return 'demux' if 'gt_demux' in step else 'align'


def demux_method(demux_type):
    """'demux_type' -> (method, add_to_previous)."""
    if demux_type is None or str(demux_type).lower() in ('', 'none', 'null'):
        return None, None
    demux_type = str(demux_type).lower()
    if demux_type == 'add_solo':
        return 'solo', 'vireo'
    if demux_type == 'add_vireo':
        return 'vireo', 'solo'
    if demux_type in ('solo', 'vireo', 'both'):
        return demux_type, None
    raise ValueError(f"Unrecognised demux_type: '{demux_type}'")


def derived(old) -> dict:
    """New path -> value, for the old keys whose meaning changed."""
    new = {}
    if 'last_step' in old:
        step = str(old['last_step'])
        new['workflow.aligner'] = aligner(step)
        new['workflow.until'] = until(step)
        new['experiment.demultiplexing.genotype.multi_vcf'] = (
            step.lower().endswith('multi_vcf'))
    if 'demux_type' in old:
        method, previous = demux_method(old['demux_type'])
        new['experiment.demultiplexing.method'] = method
        new['experiment.demultiplexing.add_to_previous'] = previous
    if 'demux_run_type' in old:
        # Only 'wo_gt' meant something: no donor genotypes (donor1, donor2, ...)
        new['experiment.demultiplexing.genotype.enabled'] = (
            str(old['demux_run_type']).lower() != 'wo_gt')
    if 'hto_demux_type' in old:
        hto = old['hto_demux_type']
        new['experiment.hto.enabled'] = hto is not None
        new['experiment.hto.multiplexing'] = (
            str(hto).lower() if hto is not None else None)
    return new


def eol_comment(cm, key):
    """The end-of-line comment of 'key' in 'cm' (not the comment lines below it)."""
    tok = cm.ca.items.get(key, [None] * 4)[2]
    if tok is None or tok.value.lstrip(' ').startswith('\n'):
        return None
    line = tok.value.split('\n', 1)[0].strip()
    return line[1:].strip() if line.startswith('#') else None


def drop_trailing_comments(node):
    """Drop the comment lines after the last key of a moved block: they
    were the header of the next block in the old layout."""
    while isinstance(node, (CommentedMap, CommentedSeq)) and len(node):
        last = list(node.keys())[-1] if isinstance(node, CommentedMap) else len(node) - 1
        items = node.ca.items.get(last)
        if items:
            pos = 2 if isinstance(node, CommentedMap) else 0
            if items[pos] is not None:
                text = (eol_comment(node, last) if pos == 2 else None)
                items[pos] = None
                if text:
                    node.yaml_add_eol_comment(text, last)
        node = node[last]


def comment_tokens(node):
    """All comment tokens in a block (its own and its children's)."""
    if not isinstance(node, (CommentedMap, CommentedSeq)):
        return
    for items in node.ca.items.values():
        yield from (t for t in items if t is not None and hasattr(t, 'value'))
    for child in (node.values() if isinstance(node, CommentedMap) else node):
        yield from comment_tokens(child)


def indent_comments(node, shift: int):
    """Indent the comment lines of a block moved 'shift' columns deeper
    (ruamel keeps their old column, the keys get the new one)."""
    for tok in comment_tokens(node):
        first, *rest = tok.value.split('\n')
        rest = [' ' * shift + l if l.strip().startswith('#') else l for l in rest]
        tok.value = '\n'.join([first] + rest)


def keep_anchors(node):
    """Keep the anchors no alias uses (e.g. other genomes to pick)."""
    if isinstance(node, (CommentedMap, CommentedSeq)):
        if node.anchor.value:
            node.anchor.always_dump = True
        for child in (node.values() if isinstance(node, CommentedMap) else node):
            keep_anchors(child)
    elif getattr(getattr(node, 'anchor', None), 'value', None):
        node.anchor.always_dump = True


def reorder(node, order):
    """Put the keys of 'order' first, in that order."""
    for key in reversed(order):
        if key in node:
            comment = node.ca.items.get(key)
            node.insert(0, key, node.pop(key))
            if comment:
                node.ca.items[key] = comment


def set_path(root, path: str, value, comment=None, first=False):
    *parents, key = path.split('.')
    node = root
    for p in parents:
        if p not in node:
            node[p] = CommentedMap()
        node = node[p]
    if key in node:
        raise ValueError(f"'{path}' is set twice")
    if first:
        # After the folder structures already there
        pos = next((i for i, k in enumerate(node) if not k.startswith('fold_struct')),
                   len(node))
        node.insert(pos, key, value)
    else:
        node[key] = value
    if comment:
        node.yaml_add_eol_comment(comment, key)


def migrate(old):
    """The new layout of 'old' (a CommentedMap of the old layout)."""
    if not isinstance(old, CommentedMap):
        raise ValueError("The config isn't a mapping")
    if any(s in old for s in ['workflow', 'experiment', 'modules', 'references']):
        raise ValueError("Already in the new layout")
    unknown = set(old) - set(MOVES) - UNUSED - {
        'last_step', 'demux_type', 'demux_run_type', 'hto_demux_type'}
    if unknown:
        raise ValueError(f"Unrecognised keys: {', '.join(sorted(unknown))}")

    keep_anchors(old)
    new = CommentedMap()
    # Blocks first, so the folder structures can be added to them
    keys = sorted(old, key=lambda k: k.startswith('fold_struct'))
    for key in keys:
        if key not in MOVES:
            continue
        value = old[key]
        drop_trailing_comments(value)
        # A block's keys were 2 columns in, now 2 per level of its new path
        indent_comments(value, 2 * (MOVES[key].count('.')))
        set_path(new, MOVES[key], value, eol_comment(old, key),
                 first=key.startswith('fold_struct'))
    for path, value in derived(old).items():
        set_path(new, path, value)

    for path, order in ORDER.items():
        node = new
        for p in path.split('.'):
            node = node.get(p) if isinstance(node, CommentedMap) else None
        if isinstance(node, CommentedMap):
            reorder(node, order)

    # Sections in a fixed order, each with a header
    out = CommentedMap()
    for s, header in SECTIONS.items():
        if s in new:
            out[s] = new[s]
            out.yaml_set_comment_before_after_key(
                s, before=('\n' if len(out) > 1 else '') + header)
    return out


def changes(old) -> list:
    """What 'migrate' does to 'old', one line each."""
    lines = [f"{k} -> {MOVES[k]}" for k in old if k in MOVES]
    lines += [f"{k}: dropped (unused)" for k in old if k in UNUSED]
    if 'demux_run_type' in old and str(old['demux_run_type']).lower() not in ('wo_gt', 'w_gt'):
        lines.append(
            f"demux_run_type '{old['demux_run_type']}': only its genotype "
            "use is kept; set run_tag if it named the run's sub-folder")
    lines += [f"{p} = {v}" for p, v in derived(old).items()]
    return lines


def yaml_rt() -> YAML:
    y = YAML()
    y.preserve_quotes = True
    y.width = 4096
    y.indent(mapping=2, sequence=4, offset=2)
    # 'null', not an empty value
    y.representer.add_representer(
        type(None),
        lambda r, _: r.represent_scalar('tag:yaml.org,2002:null', 'null'))
    return y


def migrate_text(text: str, quiet: bool = True) -> str:
    y = yaml_rt()
    try:
        old = y.load(text)
    except DuplicateKeyError as e:
        # PyYAML (snakemake) keeps the last value, ruamel the first
        raise ValueError(
            "Duplicate key: snakemake uses its last value; remove the "
            f"other one first.\n{e.problem}\n{e.problem_mark}")
    if not quiet:
        print('\n'.join(changes(old)), file=sys.stderr)
    # The comment lines the file starts with (e.g. a dry-run scenario's checks)
    header = []
    for line in text.splitlines(keepends=True):
        if line.strip() and not line.lstrip().startswith('#'):
            break
        header.append(line)
    out = io.StringIO()
    y.dump(migrate(old), out)
    return ''.join(header) + out.getvalue()


def main():
    parser = argparse.ArgumentParser(
        description="Convert configs from the old layout to the new one.")
    parser.add_argument('configs', nargs='+', help="Configs in the old layout.")
    out = parser.add_mutually_exclusive_group(required=True)
    out.add_argument('-o', '--output', help="Output (one config only).")
    out.add_argument('--in-place', action='store_true',
                     help="Overwrite each config.")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Don't list the changes.")
    args = parser.parse_args()

    if args.output and len(args.configs) > 1:
        parser.error("-o takes one config; use --in-place for several")

    failed = False
    for path in args.configs:
        if not args.quiet:
            print(f"== {path}", file=sys.stderr)
        try:
            with open(path) as f:
                text = migrate_text(f.read(), args.quiet)
        except ValueError as e:
            print(f"{path}: {e}", file=sys.stderr)
            failed = True
            continue
        with open(args.output or path, 'w') as f:
            f.write(text)
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
