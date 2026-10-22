#!/usr/bin/env python3
"""Where the demultiplexing rules write their outputs.

- Rule A (``demux_run_<tool>``, one per tool of this run): an h5ad with
  the demultiplexing of one tool (``solo``, ``vireo``, or ``none`` for no
  demultiplexing; QC stats in ``uns['qc_stats']``) and the stats of that
  tool (``_hashsolo_stats.tsv`` / ``_vireo_stats.tsv``), in its dirs.
- Rule B (``demux_final``): the final h5ad of a pool, built from the
  Rule A outputs picked for it (the swap correction file, else this run),
  its basic stats (QC and the demultiplex versions used) and its
  demultiplex info (basic and tool stats together, read by update_logs).
"""

import os
from dataclasses import dataclass, field
from functools import lru_cache
from typing import List, Optional
from .utils import get_demux_type, read_files_ext, run_dir


# Tools run by Rule A for each demux_type ('none': only create the h5ad)
TOOLS_BY_DEMUX_TYPE = {
    None: ['none'],
    'solo': ['solo'],
    'add_solo': ['solo'],
    'vireo': ['vireo'],
    'add_vireo': ['vireo'],
    'both': ['solo', 'vireo'],
}

# The tool whose outputs (from 'previous_run_tag') an add_* run adds to
PREVIOUS_TOOL = {'add_solo': 'vireo', 'add_vireo': 'solo'}


def demux_tools(config) -> List[str]:
    return TOOLS_BY_DEMUX_TYPE[get_demux_type(config)]


def is_multi_vcf(config) -> bool:
    return config['last_step'].lower().endswith('multi_vcf')


def is_multiome(config) -> bool:
    return 'multiome' in config['last_step'].lower()


def vcf_types(config) -> List[str]:
    """Values of the 'vcf_type' wildcard in multi-VCF runs: the vcf columns
    of 'vcf_info' ('gt_demux_pipeline: vcf_info_columns: vcf')."""
    if not is_multi_vcf(config):
        return []
    if '{vcf_type}' not in config['fold_struct_gt_demux']:
        raise ValueError(
            "Multi-VCF runs need '{vcf_type}' in 'fold_struct_gt_demux' "
            "(e.g. \"Sample_{pool}-cDNA/{vcf_type}/\")"
        )
    cols = config['gt_demux_pipeline']['vcf_info_columns']['vcf']
    if not cols:
        raise ValueError(
            "Multi-VCF runs need the vcf columns in "
            "'gt_demux_pipeline: vcf_info_columns: vcf'"
        )
    return [cols] if isinstance(cols, str) else list(cols)


def check_supported(config) -> None:
    """Raise for combinations of last_step and demux_type not supported yet."""
    tools = demux_tools(config)
    if is_multiome(config) and 'solo' in tools:
        raise ValueError(
            "hashsolo (demux_type solo, add_solo, both or add_vireo) isn't "
            "supported for multiome runs yet"
        )
    if config['last_step'].lower().startswith('cellranger') and 'solo' in tools:
        raise ValueError(
            "hashsolo isn't supported with cellranger count yet (its "
            "barcodes end with '-1', the HTO counts' don't)"
        )
    if is_multi_vcf(config) and tools != ['vireo']:
        raise ValueError(
            "Multi-VCF runs support only demux_type vireo for now"
        )


def pool_fold(fold: str) -> str:
    """'fold' without its '{modality}' part, so one file covers a pool."""
    parts = [p for p in fold.split('/') if p and '{modality}' not in p]
    return '/'.join(parts)


def _vcf_fold(config, fold: str) -> str:
    # Multi-VCF inputs carry the 'vcf_type' wildcard, so the outputs must too
    if is_multi_vcf(config) and '{vcf_type}' not in fold:
        return fold + "_{vcf_type}"
    return fold


# Stats file of each tool (Rule A)
STATS_SUFFIX = {'solo': '_hashsolo_stats.tsv', 'vireo': '_vireo_stats.tsv'}


@dataclass
class DemuxRun:
    """Rule A outputs of one tool in the run tagged 'run_tag'."""
    tool: str
    run_tag: Optional[str]
    h5ad: str
    stats: Optional[str] # None for 'none'

    @property
    def label(self) -> str:
        return f"{self.tool}:{self.run_tag or ''}"


def demux_run(config, tool: str, run_tag: Optional[str]) -> DemuxRun:
    if tool == 'solo':
        section, fold = 'hashsolo_demux_pipeline', config['fold_struct_demux']
    elif tool == 'vireo':
        section, fold = 'gt_demux_pipeline', pool_fold(config['fold_struct_gt_demux'])
    elif tool == 'none':
        section, fold = 'gt_demux_pipeline', config['fold_struct']
    else:
        raise ValueError(f"Unknown demultiplexing tool: '{tool}'")

    fold = _vcf_fold(config, fold)
    h5ad = config[section]['final_count_matrix_h5ad']
    if tool == 'none':
        h5ad = h5ad.replace("_vS", "")

    stats = None
    if tool in STATS_SUFFIX:
        stats = (
            f"{run_dir(config, section, 'demultiplex_info_dir', run_tag)}"
            f"{fold}{STATS_SUFFIX[tool]}"
        )

    return DemuxRun(
        tool=tool,
        run_tag=run_tag,
        h5ad=f"{run_dir(config, section, 'final_count_matrix_dir', run_tag)}{fold}{h5ad}",
        stats=stats,
    )


def parse_label(label: str) -> List[tuple]:
    """'vireo:w_gt_2' -> [('vireo', 'w_gt_2')]; 'both:x' -> solo and vireo.

    No run_tag ('vireo' or 'vireo:') means the untagged dirs. Runs of
    different tools are separated by commas: 'solo:hto,vireo:w_gt'.
    """
    if ',' in str(label):
        return [run for part in str(label).split(',') for run in parse_label(part)]
    tool, _, run_tag = str(label).strip().partition(':')
    tool = tool.strip().lower()
    run_tag = run_tag.strip() or None
    tools = ['solo', 'vireo'] if tool == 'both' else [tool]
    for t in tools:
        if t not in ('solo', 'vireo', 'none'):
            raise ValueError(
                f"Unknown tool '{t}' in demultiplex version '{label}' "
                "(expected solo, vireo or both, then ':<run_tag>')"
            )
    return [(t, run_tag) for t in tools]


@lru_cache(maxsize=None)
def _read_swap_file(fname: str):
    return read_files_ext(fname)


def swap_correction(config) -> dict:
    return config['demultiplex'].get('swap_correction') or {}


def swap_version(config, pool: str) -> Optional[str]:
    """Demultiplex version of 'pool' in the swap correction file, if any."""
    swap = swap_correction(config)
    if not swap.get('file'):
        return None
    if not os.path.isfile(swap['file']):
        raise ValueError(
            "'demultiplex: swap_correction: file' doesn't exist: "
            f"{swap['file']}"
        )
    df = _read_swap_file(swap['file'])
    versions = df.loc[
        df[swap['pool_column']].astype(str) == str(pool),
        swap['version_column']
    ].dropna().unique()
    if len(versions) > 1:
        raise ValueError(
            f"Pool '{pool}' has more than one demultiplex version in "
            f"{swap['file']}: {list(versions)}"
        )
    return versions[0] if len(versions) else None


def final_sources(config, pool: str) -> List[DemuxRun]:
    """Rule A outputs that make up the final h5ad of 'pool', oldest first."""
    version = swap_version(config, pool)
    if version is not None:
        return [demux_run(config, t, r) for t, r in parse_label(version)]

    demux_type = get_demux_type(config)
    ret = []
    if demux_type in PREVIOUS_TOOL:
        ret.append(demux_run(
            config, PREVIOUS_TOOL[demux_type],
            config['demultiplex'].get('previous_run_tag')
        ))
    ret.extend(
        demux_run(config, t, config.get('run_tag')) for t in demux_tools(config)
    )
    return ret


@dataclass
class FinalPaths:
    """Rule B outputs, and the table of versions used for all pools."""
    h5ad: str
    basic_stats: str
    info: str
    versions: str


def final_paths(config) -> FinalPaths:
    c = config['demultiplex']
    if not c.get('demux_count_matrix_dir'):
        raise ValueError(
            "Set 'demultiplex: demux_count_matrix_dir', where the final "
            "h5ad of each pool is written"
        )
    fold = _vcf_fold(config, config['fold_struct_final'])
    return FinalPaths(
        h5ad=f"{c['demux_count_matrix_dir']}{fold}{c['final_count_matrix_h5ad']}",
        basic_stats=f"{c['demultiplex_info_dir']}{fold}_basic_stats.tsv",
        info=f"{c['demultiplex_info_dir']}{fold}{c['demultiplex_info']}",
        versions=(
            f"{c['demultiplex_info_dir']}final_demux_versions"
            f"{'_{vcf_type}' if is_multi_vcf(config) else ''}.tsv"
        ),
    )


def raw_output(config, tool: str) -> str:
    """Output of the demultiplexing tool itself (before Rule A)."""
    if tool == 'vireo':
        c = config['gt_demux_pipeline']
        return f"{c['vireosnp_dir']}{config['fold_struct_gt_demux']}{c['donors_classification']}"
    c = config['hashsolo_demux_pipeline']
    return f"{c['calico_solo_dir']}{config['fold_struct_demux']}{c['calico_solo_h5ad']}"


@dataclass
class DemuxConfig:
    """Demultiplexing targets of this run."""

    targets: List[str] = field(default_factory=list)
    # Not per pool (see 'collect_demux_versions')
    summary_targets: List[str] = field(default_factory=list)
    multiome: bool = False
    create_h5ad: bool = False

    @classmethod
    def from_snakemake(cls, config):
        check_supported(config)
        create_h5ad = bool(config.get("create_final_h5ad", False))
        tools = demux_tools(config)

        summary_targets = []
        if create_h5ad or tools == ['none']:
            targets = [final_paths(config).h5ad]
            summary_targets = [final_paths(config).versions]
        else:
            # Stop at the outputs of the demultiplexing tools
            targets = [raw_output(config, t) for t in tools]

        return cls(
            targets=targets,
            summary_targets=summary_targets,
            multiome=is_multiome(config),
            create_h5ad=create_h5ad,
        )
