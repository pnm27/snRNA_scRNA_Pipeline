#!/usr/bin/env python3
"""Convert a config from the old (flat) layout to the new one.

New layout::

  references  genomes, barcodes, regions and tool versions
  workflow    aligner, start, until (replace 'last_step'), picard_metrics,
              gt_check
  experiment  inputs, hto, demultiplexing (replaces 'demux_type',
              'demux_run_type' and 'hto_demux_type'), qc
  modules     the '*_pipeline' blocks, each with its own 'fold_struct*'
  execution   rule resources (was 'resources')

New keys are written out with their values for the old configs:
'workflow: start: align' (wherever 'last_step' is converted) and, in a
whole config (one with a whole split bams block), 'split_by: external'
(donor labels from outside the pipeline, unused unless 'split_by: input:
external'), the WASP block 'modules: wasp' and its env module
'references: tools: wasp' (unused unless 'until: wasp'). A last_step
naming wasp (e.g. STARsolo_split_bams_wasp) becomes 'until: wasp'.

A config already in the new layout gets only the keys added since that
it lacks (e.g. 'modules: wasp', 'split_bams: modalities'), and the WASP
block's 'remap' and 'paired_end' move under 'rna' (with ATAC WASP); it's
left as is when it has them all.

Comments inside the moved blocks and the anchors (e.g. ``*3v3``) are kept.
Configs that only override some keys (e.g. the dry-run scenarios) are
converted too: only the keys they have are moved.

Usage::

    python3 migrate_config.py new_config.yaml -o new_config.migrated.yaml
    python3 migrate_config.py tests/dryrun/scenarios/*.yaml --in-place
"""

import argparse, io, sys
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString
from ruamel.yaml.constructor import DuplicateKeyError
from ruamel.yaml.error import CommentMark
from ruamel.yaml.tokens import CommentToken


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
    'workflow':   "What to run: the aligner, the first stage ('start') and the last ('until')",
    'experiment': "What was sequenced and how it's multiplexed",
    'modules':    "Parameters, output dirs and folder structures of each module",
    'execution':  "Rule resources (override rules/resources.smk)",
}

# Order of the keys made by the migration (other keys follow)
ORDER = {
    'workflow':                 ['aligner', 'start', 'until', 'picard_metrics', 'gt_check'],
    'experiment':               ['inputs', 'hto', 'demultiplexing', 'qc'],
    'experiment.demultiplexing': ['method', 'add_to_previous', 'run_tag',
                                  'genotype', 'create_final_h5ad'],
    'experiment.demultiplexing.genotype': ['enabled', 'multi_vcf'],
    'experiment.hto':           ['enabled', 'multiplexing'],
}

# The stages 'until' can name, in order
STAGES = ['align', 'cellsnp', 'demux', 'split_bams', 'identify_swaps', 'wasp']

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
    for stage in ['wasp', 'identify_swaps', 'split_bams']:
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


def external_labels_block() -> CommentedMap:
    """'split_by: external' with its defaults (see new_config.example.yaml)."""
    def step(**kw):
        m = CommentedMap(type='pythonic', **kw)
        m.fa.set_flow_style()
        return m
    ext = CommentedMap()
    ext['file'] = None
    ext['donor_column'] = 'individualID'
    ext['pool_column'] = None
    ext['barcode_column'] = None
    ext['barcode_transformations'] = CommentedSeq([
        step(operation='split', by='_'),
        step(operation='select_index', index=-1),
    ])
    pools = CommentedSeq()
    pools.fa.set_flow_style()
    ext['pool_transformations'] = pools
    ext['pools'] = None
    ext['min_barcode_fraction'] = 0.02
    ext.yaml_add_eol_comment('pooled h5ad, or csv/tsv/txt(.gz) table, one row per cell', 'file')
    ext.yaml_add_eol_comment('obs (h5ad) or table column with the pool of origin', 'pool_column')
    ext.yaml_add_eol_comment('null: obs names (h5ad only)', 'barcode_column')
    ext.yaml_add_eol_comment(
        "null: the aligner's bams_dir; or a table (pool, bam[, bai])", 'pools')
    ext.yaml_add_eol_comment(
        "least fraction of sampled CB reads from the labelled cells", 'min_barcode_fraction')
    ext.yaml_set_comment_before_after_key(
        'barcode_transformations',
        before="Steps from the file's unique barcodes to the bam's CB tags: the "
               "part after the last '_'", indent=8)
    return ext


# Env module of WASP written to 'references: tools: wasp' (as new_config.example.yaml)
WASP_MODULE = 'wasp/2021-04'

SAMPLE_NAMES_COMMENT = ("pool[, donor], final_donor, vcfID: final donors and genotype ID(s) "
                        "of pools (or donors)")


RNA_COMMENT = "cDNA bams: STAR remapping, duplicates kept"
ATAC_COMMENT = "ATAC bams (multiome): paired-end, bwa-mem remapping"
MODALITIES_COMMENT = ("bams split per donor; multiome: [cDNA, ATAC] (both by the "
                      "cDNA or external labels)")


def eol(text: str) -> CommentToken:
    """An end-of-line comment token, two spaces after the value (as in the
    example config)."""
    return CommentToken(f" # {text}\n", CommentMark(0), None)


def set_eol(node, key, text):
    node.ca.items[key] = [None, None, eol(text), None]


def rna_block() -> CommentedMap:
    """'modules: wasp: rna' with its defaults."""
    remap = CommentedMap()
    remap['genome'] = None
    remap['gtf'] = None
    remap['extra_params'] = None
    remap.yaml_add_eol_comment(
        "STAR index; null: the aligner's reference (start: wasp: required)", 'genome')
    remap.yaml_add_eol_comment(
        "null: STARsolo's genome_pick gtf when the genome is the aligner's", 'gtf')
    rna = CommentedMap()
    rna['remap'] = remap
    rna['paired_end'] = False
    rna.yaml_add_eol_comment("true for paired-end GEX reads; 10x 3': false", 'paired_end')
    return rna


def atac_block() -> CommentedMap:
    """'modules: wasp: atac' with its defaults."""
    remap = CommentedMap()
    remap['bwa_index'] = None
    remap['extra_params'] = DoubleQuotedScalarString("-M -Y")
    set_eol(remap, 'bwa_index', "bwa-indexed fasta; null: cellranger-arc's "
            "fasta/genome.fa (start: wasp: required)")
    atac = CommentedMap()
    atac['remap'] = remap
    atac['drop_duplicates'] = 'marked'
    set_eol(atac, 'drop_duplicates', "marked: remove the aligner's duplicates "
            "(0x400) after WASP; or none")
    return atac


def modalities_seq() -> CommentedSeq:
    seq = CommentedSeq(['cDNA'])
    seq.fa.set_flow_style()
    return seq


def wasp_block() -> CommentedMap:
    """'modules: wasp' with its defaults (see new_config.example.yaml)."""
    wasp = CommentedMap()
    wasp['samples'] = None
    wasp['sample_names'] = None
    wasp['genotypes'] = None
    wasp['rna'] = rna_block()
    wasp['atac'] = atac_block()
    wasp['wasp_dir'] = 'projdir/wasp_bams/'
    wasp['work_dir'] = 'projdir/wasp_tmp/'
    wasp.yaml_set_comment_before_after_key(
        'samples',
        before="pool, donor (as split), final_donor (after swap correction; names the "
               "outputs), vcfID\n(genotype ID(s), ';'- or ','-separated; "
               "find_intersecting_snps.py's --samples), modality\n(cDNA or ATAC; empty: "
               "cDNA with start: wasp, else every split modality), and bam[, bai]\nwith "
               "start: wasp; genotype_set with several sets. Empty or no final_donor / "
               "vcfID: from\nsample_names, else the donor (no final_donor column in "
               "either table: a warning)", indent=4)
    wasp.yaml_add_eol_comment(SAMPLE_NAMES_COMMENT, 'sample_names')
    wasp.yaml_add_eol_comment(
        "sets, e.g. main: {h5_dir: path/} or main: {vcf: x.bcf, chrom_sizes: chrom.txt}",
        'genotypes')
    set_eol(wasp, 'rna', RNA_COMMENT)
    set_eol(wasp, 'atac', ATAC_COMMENT)
    wasp.yaml_add_eol_comment(
        'large temporary files, removed per job', 'work_dir')
    return wasp


def whole_config(old) -> bool:
    """A whole config (not one overriding some keys): it has a whole split
    bams block."""
    sb = old.get('split_bams_pipeline')
    return isinstance(sb, dict) and 'inp_split_bams_dir' in sb


def new_layout(cfg) -> bool:
    """In the new layout: it has a new section and no old key."""
    sections = [k for k in cfg if k in SECTIONS]
    if not sections:
        return False
    old = [k for k in cfg if k not in SECTIONS]
    if old:
        raise ValueError(
            "Mixes the old and new layouts: new sections "
            f"({', '.join(sections)}) and old keys ({', '.join(old)})")
    return True


def added_since(cfg) -> dict:
    """New path -> value, for the keys added since the new layout that a
    config in the new layout lacks (as 'derived' writes them out for old
    configs)."""
    def get(path):
        node = cfg
        for p in path.split('.'):
            if not isinstance(node, dict) or p not in node:
                return None, False
            node = node[p]
        return node, True

    new = {}
    sb, _ = get('modules.split_bams')
    whole = isinstance(sb, dict) and 'inp_split_bams_dir' in sb
    if whole:
        split_by = sb.get('split_by')
        if isinstance(split_by, dict) and 'external' not in split_by:
            new['modules.split_bams.split_by.external'] = external_labels_block()
        if not get('references.tools.wasp')[1]:
            new['references.tools.wasp'] = WASP_MODULE
        if not get('modules.wasp')[1]:
            new['modules.wasp'] = wasp_block()
    ext, _ = get('modules.split_bams.split_by.external')
    if isinstance(ext, dict) and 'min_barcode_fraction' not in ext:
        new['modules.split_bams.split_by.external.min_barcode_fraction'] = 0.02
    if whole and 'modalities' not in sb:
        new['modules.split_bams.modalities'] = modalities_seq()
    wasp, _ = get('modules.wasp')
    if isinstance(wasp, dict) and 'samples' in wasp and 'sample_names' not in wasp:
        new['modules.wasp.sample_names'] = None
    # 'remap' and 'paired_end' (cDNA only, before ATAC WASP) go to 'rna'
    if isinstance(wasp, dict) and 'rna' not in wasp and (
            'remap' in wasp or 'paired_end' in wasp):
        new['modules.wasp.rna'] = MOVE_TO_RNA
    if isinstance(wasp, dict) and 'wasp_dir' in wasp and 'atac' not in wasp:
        new['modules.wasp.atac'] = atac_block()
    swap, _ = get('modules.demultiplex.swap_correction')
    if isinstance(swap, dict) and 'donor_to' in swap:
        for key in ('status_column', 'versions'):
            if key not in swap:
                new[f'modules.demultiplex.swap_correction.{key}'] = None
    return new


def last_leaf(node):
    """(map, key) of a block's deepest last key, which holds the comment
    lines after the block."""
    while isinstance(node, CommentedMap) and len(node):
        key = list(node)[-1]
        if not (isinstance(node[key], CommentedMap) and len(node[key])):
            return node, key
        node = node[key]
    return None, None


def take_tail(node) -> str:
    """Remove the comment lines after a block (e.g. the next section's
    header), keeping its last key's end-of-line comment; return them."""
    m, k = last_leaf(node)
    items = m.ca.items.get(k) if m is not None else None
    if not items or items[2] is None:
        return ''
    first, _, rest = items[2].value.partition('\n')
    if first.strip().startswith('#'):
        items[2].value = first + '\n'
    else:
        items[2] = None
    return rest


def put_tail(node, tail: str):
    """Put comment lines after a block (see 'take_tail')."""
    m, k = last_leaf(node)
    if not tail or m is None:
        return
    items = m.ca.items.setdefault(k, [None] * 4)
    if items[2] is not None:
        items[2].value = items[2].value.rstrip('\n') + '\n' + tail
    else:
        items[2] = CommentToken('\n' + tail, CommentMark(0), None)


# Keys of 'modules: wasp' moved under 'rna' (see added_since)
MOVE_TO_RNA = ('remap', 'paired_end')


def move_to_rna(wasp):
    """Move 'remap' and 'paired_end' of a WASP block under 'rna', where
    they were; their comments too."""
    keys = [k for k in MOVE_TO_RNA if k in wasp]
    pos = list(wasp).index(keys[0])
    rna = CommentedMap()
    for k in keys:
        comment = wasp.ca.items.pop(k, None)
        rna[k] = wasp.pop(k)
        if comment:
            rna.ca.items[k] = comment
    indent_comments(rna, 2)
    for tok in comment_tokens(rna):
        # end-of-line comments: two columns further, as their keys
        if not tok.value.startswith('\n'):
            tok.start_mark.column += 2
    wasp.insert(pos, 'rna', rna)
    set_eol(wasp, 'rna', RNA_COMMENT)


def update(cfg) -> dict:
    """Add the keys of 'added_since' to 'cfg' (in the new layout), in
    place; return them."""
    new = added_since(cfg)
    for path, value in new.items():
        parent = cfg
        for p in path.split('.')[:-1]:
            parent = parent[p]
        if path == 'modules.wasp.rna':
            move_to_rna(parent)
            continue
        if path == 'modules.wasp.atac':
            # After 'rna' (or at the end)
            pos = list(parent).index('rna') + 1 if 'rna' in parent else len(parent)
            parent.insert(pos, 'atac', value)
            set_eol(parent, 'atac', ATAC_COMMENT)
            continue
        if path == 'modules.split_bams.modalities':
            pos = (list(parent).index('min_cells_per_donor')
                   if 'min_cells_per_donor' in parent else len(parent))
            parent.insert(pos, 'modalities', value)
            set_eol(parent, 'modalities', MODALITIES_COMMENT)
            continue
        if path == 'modules.wasp.sample_names':
            # Next to the samples table
            parent.insert(list(parent).index('samples') + 1, 'sample_names', value)
            # (two spaces before the '#', as in the example config)
            parent.ca.items['sample_names'] = [None, None, CommentToken(
                f" # {SAMPLE_NAMES_COMMENT}\n", CommentMark(0), None), None]
            continue
        tail = take_tail(parent)
        set_path(cfg, path, value)
        put_tail(parent, tail)
    return new


def derived(old) -> dict:
    """New path -> value, for the old keys whose meaning changed, and the
    new keys written out for old configs."""
    new = {}
    if 'last_step' in old:
        step = str(old['last_step'])
        new['workflow.aligner'] = aligner(step)
        new['workflow.start'] = 'align'  # from the fastqs, as before
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
    # A whole split bams block (not a config overriding some of its keys)
    sb = old.get('split_bams_pipeline')
    if (whole_config(old) and isinstance(sb.get('split_by'), dict)
            and 'external' not in sb['split_by']):
        new['modules.split_bams.split_by.external'] = external_labels_block()
    # The WASP stage, added since (unused unless 'until: wasp')
    if whole_config(old):
        new['references.tools.wasp'] = WASP_MODULE
        new['modules.wasp'] = wasp_block()
        new['modules.split_bams.modalities'] = modalities_seq()
    # Keys added to a whole swap correction block since (see
    # new_config.example.yaml)
    swap = (old.get('demultiplex') or {}).get('swap_correction')
    if isinstance(swap, dict) and 'donor_to' in swap:
        for key in ('status_column', 'versions'):
            if key not in swap:
                new[f'modules.demultiplex.swap_correction.{key}'] = None
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


def get_node(root, path: str):
    for p in path.split('.'):
        root = root[p]
    return root


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
    if new_layout(old):
        raise ValueError("Already in the new layout (see 'update')")
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
        if path == 'modules.split_bams.modalities':
            set_eol(get_node(new, path.rsplit('.', 1)[0]), 'modalities', MODALITIES_COMMENT)

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
    lines += [f"{p} = {v}" if not isinstance(v, dict) else f"{p}: added (defaults)"
              for p, v in derived(old).items()]
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
    if isinstance(old, CommentedMap) and new_layout(old):
        # Already converted: only add the keys added since
        keep_anchors(old)
        added = update(old)
        if not quiet:
            print('\n'.join(f"{p}: added" for p in added)
                  or "In the new layout, nothing to add", file=sys.stderr)
        if not added:
            return text
        out = io.StringIO()
        y.dump(old, out)
        return out.getvalue()
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
        description="Convert configs from the old layout to the new one, "
                    "or add the keys added since to configs in the new layout.")
    parser.add_argument('configs', nargs='+', help="Configs (old or new layout).")
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
