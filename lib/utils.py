#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Any, Optional, Union # Need verion > 3.10
import pandas as pd, os


@dataclass(frozen=True)
class RunSettings:
    """What this run does: the 'workflow' and 'experiment' keys the rules
    branch on (checked by schema/config.schema.yaml), read in one place."""
    aligner: Optional[str]          # starsolo, cellranger or cellranger_arc
    until: Optional[str]            # last stage (rules/produce_targets.smk)
    method: Optional[str]           # solo, vireo, both or None
    add_to_previous: Optional[str]  # tool of the earlier run added to
    run_tag: Optional[str]
    genotyped: bool                 # else vireo names donors donor0, donor1, ...
    multi_vcf: bool
    hto_multiplexing: Optional[str] # single or multiple; None without HTOs

    @property
    def demux_type(self) -> Optional[str]:
        """solo, vireo, both, add_solo / add_vireo (method added to a run of
        the other tool), or None when no demux is needed."""
        if self.method is None:
            return None
        return f"add_{self.method}" if self.add_to_previous else self.method


def _lower(value) -> Optional[str]:
    if value is None or str(value).lower() in ('none', 'null', ''):
        return None
    return str(value).lower()


def run_settings(config) -> RunSettings:
    wf = config.get('workflow') or {}
    exp = config.get('experiment') or {}
    d = exp.get('demultiplexing') or {}
    gt = d.get('genotype') or {}
    hto = exp.get('hto') or {}
    return RunSettings(
        aligner=_lower(wf.get('aligner')),
        until=_lower(wf.get('until')),
        method=_lower(d.get('method')),
        add_to_previous=_lower(d.get('add_to_previous')),
        run_tag=d.get('run_tag'),
        genotyped=bool(gt.get('enabled')),
        multi_vcf=bool(gt.get('multi_vcf')),
        hto_multiplexing=_lower(hto.get('multiplexing')) if hto.get('enabled') else None,
    )


def get_demux_type(config) -> Union[None, str]:
    """See RunSettings.demux_type."""
    return run_settings(config).demux_type


def is_genotyped(config) -> bool:
    """Whether donor genotypes are used (else vireo names donors donor0, donor1, ...)."""
    return run_settings(config).genotyped


# Output dirs of a demultiplexing run; 'run_tag' adds a sub-folder to each
RUN_TAGGED_DIRS = {
    'hashsolo': [
        'calico_solo_dir', 'final_count_matrix_dir', 'demultiplex_info_dir',
    ],
    'gt_demux': [
        'inp_for_cellsnp_dir', 'cellsnp_dir', 'vireosnp_dir',
        'final_count_matrix_dir', 'demultiplex_info_dir',
    ],
}


def apply_run_tag(config) -> None:
    """Add the 'run_tag' sub-folder to the dirs in RUN_TAGGED_DIRS.

    The untagged dirs are kept in config['_untagged_dirs'], so the outputs
    of other runs can still be found (see run_dir). Calling it again
    doesn't add the tag twice.
    """
    if '_untagged_dirs' in config:
        return

    config['_untagged_dirs'] = {
        section: {key: config['modules'][section][key] for key in keys}
        for section, keys in RUN_TAGGED_DIRS.items()
    }
    for section, keys in RUN_TAGGED_DIRS.items():
        for key in keys:
            config['modules'][section][key] = run_dir(
                config, section, key, run_settings(config).run_tag
            )


# Dirs of the split bams pipeline; bams split for genotype checks
# (gt_check) and the final per-donor bams each get a sub-folder
SPLIT_PURPOSE_DIRS = [
    'inp_split_bams_dir', 'split_bams_dir', 'per_donor_split_log_dir',
    'new_temp_dir',
]


def split_purpose(config) -> str:
    return 'gt_check' if config['workflow']['gt_check'] else 'final'


def apply_split_purpose(config) -> None:
    """Add the split_purpose sub-folder to the dirs in SPLIT_PURPOSE_DIRS."""
    if '_split_purpose' in config:
        return
    config['_split_purpose'] = split_purpose(config)
    sb = config['modules']['split_bams']
    for key in SPLIT_PURPOSE_DIRS:
        sb[key] = f"{sb[key]}{config['_split_purpose']}/"


def run_dir(config, section: str, key: str, run_tag: Union[None, str]) -> str:
    """config['modules'][section][key] of the run tagged 'run_tag' (None: untagged)."""
    base = config['_untagged_dirs'][section][key]
    return f"{base}{run_tag}/" if run_tag else base


# def read_files_ext(fname, h: Union[None, list]=None) -> pd.DataFrame :
def read_files_ext(fname: str, **kwargs: Any) -> pd.DataFrame :    
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 

    if fname.endswith('.csv'):
        return pd.read_csv(fname, **kwargs)
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t', **kwargs)
    elif fname.endswith('.txt'):
        return pd.read_csv(fname, sep=' ', **kwargs)
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")


def ret_cols(fname: str) -> int :    
    if not os.path.isfile(fname):
        raise OSError(f"The given file {fname} doesn't exist and annotations are impossible without this file!") 
    if fname.endswith('.csv'):
        return pd.read_csv(fname, nrows=1, header=None).shape[1]
    elif fname.endswith('.tsv'):
        return pd.read_csv(fname, sep='\t', nrows=1, header=None).shape[1]
    else:
        raise OSError(f"The given file {fname} doen't have either csv or tsv extension. Other extensions are not supported!")
    
