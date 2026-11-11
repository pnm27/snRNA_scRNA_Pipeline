import os, sys
from lib.config_models import DemuxConfig, vcf_types, swap_correction, swap_pools
from lib.aligner import bam_paths, get_aligner, out_prefix
from lib.utils import run_settings
from lib.io import is_external, external_pools_file, missing_pool_bams, read_pools_table
from lib.aligner import count_matrix as aligner_count_matrix


# Targets of 'rule all' come from the stages run up to 'workflow: until':
#   align           alignment (Picard's metrics instead, when 'picard_metrics'
#                   is set), always run
#   cellsnp         cellSNP's pileup (only when until is cellsnp)
#   demux           demultiplexing (see DemuxConfig)
#   split_bams      per-donor bams      } known once the donors of each pool
#   identify_swaps  QTLtools-mbv        } are (checkpoint 'create_inp_splitBams')
#
# Donor labels from outside the pipeline ('split_by: input: external') need
# no demux stage, and with 'workflow: start: split_bams' the bams are
# already there (no align stage); the pools are then the ones in the labels
# (checkpoint 'normalize_external_labels'). Split by the final h5ad ('input:
# final'), such a run makes only the final h5ads of the pools of the swap
# correction file, from the demultiplexing versions it names

# Stages run for each value of 'until'
UNTIL_STAGES = {
    "align":          ["align"],
    "cellsnp":        ["align", "cellsnp"],
    "demux":          ["align", "demux"],
    "split_bams":     ["align", "demux", "split_bams"],
    "identify_swaps": ["align", "demux", "split_bams", "identify_swaps"],
}


# Stages skipped for each value of 'start'
START_SKIPS = {
    "align":      [],
    "split_bams": ["align", "cellsnp", "demux"],
}


def until_stages(until: str, start: str = "align", external: bool = False) -> list:
    """Stages run from 'workflow: start' to 'workflow: until'."""
    if until not in UNTIL_STAGES:
        raise ValueError(
            f"Unrecognised 'workflow: until': '{until}' (expected one of "
            f"{', '.join(UNTIL_STAGES)})"
        )
    if start not in START_SKIPS:
        raise ValueError(
            f"Unrecognised 'workflow: start': '{start}' (expected one of "
            f"{', '.join(START_SKIPS)})"
        )
    skip = set(START_SKIPS[start])
    if external and until in ("split_bams", "identify_swaps"):
        skip.add("demux")  # donors come from the external labels
    return [s for s in UNTIL_STAGES[until] if s not in skip]


# ---------------------------------------------------------------------------
# This run
# ---------------------------------------------------------------------------
UNTIL = run_settings(config).until
START = run_settings(config).start
EXTERNAL = is_external(config)
POOLS_TABLE = external_pools_file(config)
STAGES = until_stages(UNTIL, START, EXTERNAL)
CFG = DemuxConfig.from_snakemake(config)
MODALITY = ["cDNA", "ATAC"] if CFG.multiome else ["cDNA"]
VCF_TYPES = vcf_types(config)

if VCF_TYPES and any(s in STAGES for s in ["split_bams", "identify_swaps"]):
    raise ValueError(
        "Splitting bams (and identifying swaps) isn't supported for "
        "multi-VCF runs yet"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# STATIC, pool-only expansion (and vcf_type, for multi-VCF runs)
def _expand_pools(templates: list) -> list:
    """Expand path templates over all pools and modalities."""
    # (an empty vcf_type list would expand to nothing)
    extra = {"vcf_type": VCF_TYPES} if VCF_TYPES else {}
    return expand(templates, pool=POOLS, modality=MODALITY, **extra)


# DYNAMIC, pool+donors expansion
def _expand_pools_donors(template: str, wildcards) -> list:
    """Expand a path template over pool + its donors, zipped."""
    wc_dict = get_real_donors({
                "pool": wildcards.pool,
                "modality": wildcards.modality
                }
            )
    if getattr(wildcards, "modality", None) not in (None, ""):
        # One value per donor: zip stops at its shortest list
        wc_dict["modality"] = [wildcards.modality] * len(wc_dict["donor"])
    return expand(template, zip, **wc_dict)


def _picard_targets() -> list:
    """Picard's outputs for 'picard_metrics' (all, gc or rnaseq; null: none)."""
    progs = (config["workflow"].get("picard_metrics") or "").lower()
    if not progs or get_aligner(config) != "starsolo":
        return []
    files = {
        "rnaseq": config["modules"]["picard"]["rnaseq_metrics"],
        "gc":     config["modules"]["picard"]["gc_summary_metrics"],
    }
    keys = files.keys() if progs == "all" else [progs]
    return [f"{out_prefix(config)}{files[k]}" for k in keys]


_NOTED = set()  # print each note once (split_pools runs more than once)


def _note_once(key: str, msg: str):
    if key not in _NOTED:
        _NOTED.add(key)
        print(msg, file=sys.stderr)


def existing_bam_pools(pools: list, source: str, hint: str) -> list:
    """Of 'pools' (those of 'source'), the ones whose bams are split in a
    'workflow: start: split_bams' run: with a pools table, the ones it lists
    (a warning for the others, e.g. bams not kept for disk space). Their
    bams must be there already: nothing aligns them in this run."""
    if POOLS_TABLE:
        table = read_pools_table(POOLS_TABLE)
        unlisted = [p for p in pools if p not in table]
        pools = [p for p in pools if p in table]
        if not pools:
            raise ValueError(
                f"None of the pools of {source} are in the pools table "
                f"{POOLS_TABLE}: {', '.join(unlisted)}. Check the table's "
                f"'pool' column{hint}"
            )
        if unlisted:
            _note_once("unlisted", (
                f"Warning: pool(s) of {source} not in the pools table "
                f"{POOLS_TABLE}, not split: {', '.join(unlisted)}"))
        extra = sorted(set(table) - set(pools))
        if extra:
            _note_once("extra", (
                f"Pools table {POOLS_TABLE}: pool(s) not in {source}, not "
                f"split: {', '.join(extra)}"))
    missing = missing_pool_bams(config, pools)
    if missing:
        where = (
            f"check the pools table ('split_by: external: pools')"
            if POOLS_TABLE else
            f"expected at {bam_paths(config, 'cDNA')['bam'].format(pool='<pool>')}: "
            "check the aligner's bams_dir and fold_struct (or list the "
            "bams in 'split_by: external: pools')"
        )
        raise ValueError(
            f"No bam for pool(s) of {source}:\n  " + "\n  ".join(missing)
            + f"\n'workflow: start: split_bams' splits existing bams: {where}{hint}"
        )
    return pools


def split_pools() -> list:
    """Pools whose bams are split: those of the external labels
    (checkpoint 'normalize_external_labels'), else the run's pools."""
    if not EXTERNAL:
        return POOLS
    import pandas as pd
    labels = checkpoints.normalize_external_labels.get().output[0]
    pools = sorted(pd.read_csv(labels, sep="\t", usecols=["pool"], dtype=str)["pool"].unique())
    if START == "split_bams":
        pools = existing_bam_pools(
            pools, f"the external labels ({labels})", ", or 'pool_transformations'")
    return pools


def swap_file_pools() -> list:
    """Pools of a 'workflow: start: split_bams' run split by the final h5ad
    ('split_by: input: final'): those of the swap correction file with a
    donor kept ('status_column')."""
    swap = swap_correction(config)
    pools = swap_pools(config)
    if swap.get('status_column'):
        dropped = sorted(set(swap_pools(config, kept_only=False)) - set(pools))
        if dropped:
            _note_once("no_status", (
                f"Warning: pool(s) of the swap correction file {swap['file']} "
                f"with no donor whose '{swap['status_column']}' is True, not "
                f"split: {', '.join(dropped)}"))
    if not pools:
        raise ValueError(
            f"No pool of the swap correction file {swap['file']} has a donor to "
            "split" + (f" ('{swap['status_column']}' True)" if swap.get('status_column') else "")
        )
    return existing_bam_pools(
        pools, f"the swap correction file ({swap['file']})", "")

# Final h5ads of existing bams ('workflow: start: split_bams', 'split_by:
# input: final'): the pools are those of the swap correction file
if START == "split_bams" and not EXTERNAL:
    POOLS = swap_file_pools()


# ---------------------------------------------------------------------------
# Per-stage target builders (path templates)
# ---------------------------------------------------------------------------

def targets_align() -> list:
    if get_aligner(config) == "cellranger_arc":
        out_dir = config['modules']['cellranger_arc']['bams_dir']
        return [
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/{f}"
            for f in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]
        ]
    if get_aligner(config) == "cellranger":
        return [aligner_count_matrix(config), bam_paths(config, "cDNA")["bai"]]
    # Picard's outputs stand for the alignment when they are asked for
    return _picard_targets() or [bam_paths(config, "cDNA")["bai"]]


def targets_cellSNP() -> list:
    c = config["modules"]["gt_demux"]
    return [f"{c['cellsnp_dir']}{config['modules']['gt_demux']['fold_struct']}{c['cellsnp_cells']}"]


def targets_demux() -> list:
    """Demultiplexing targets of this run (see DemuxConfig): the final h5ad
    of each pool (and the table of demultiplex versions used) or, without
    'create_final_h5ad', the demultiplexing tools' own outputs."""
    return list(CFG.targets) + list(CFG.summary_targets)


def targets_SplitBams() -> list:
    """Per-donor bams (split_bams_dir has the gt_check/ or final/ sub-folder)."""
    out_dir  = f"{config['modules']['split_bams']['split_bams_dir']}"
    op = config["modules"]["split_bams"]["fold_struct_bams"]
    return [
        os.path.join(out_dir, op, "{donor}.bam")
    ]


def targets_identify_swaps() -> list:
    out_dir = config["modules"]["identify_swaps"]["mbv_out_dir"]
    fs      = config["modules"]["identify_swaps"]["fold_struct"]
    suff    = config["modules"]["identify_swaps"]["mbv_suffix"]
    return [os.path.join(out_dir, f"{fs}{suff}")]


STATIC_TARGETS = {
    "align":   targets_align,
    "cellsnp": targets_cellSNP,
    "demux":   targets_demux,
}


# ---------------------------------------------------------------------------
# produce_targets  — static expansion (DAG construction time)
# ---------------------------------------------------------------------------

def produce_targets(wildcards) -> list:
    targets = []
    for stage in STAGES:
        if stage in STATIC_TARGETS:
            targets += _expand_pools(STATIC_TARGETS[stage]())

    # Per-donor targets are resolved by 'resolve_pool_targets' (Snakefile)
    if "split_bams" in STAGES and EXTERNAL:
        # Pools known once the labels are read ('resolve_external_pools'),
        # so the other targets don't wait for them
        targets.append("resolved/external_pools.done")
    elif "split_bams" in STAGES:
        targets += expand("resolved/{pool}/{modality}.done", pool=POOLS, modality=MODALITY)

    return targets


# ---------------------------------------------------------------------------
# produce_targets_dynamic — runtime expansion (after hash file exists)
# ---------------------------------------------------------------------------

def external_pool_targets(wildcards) -> list:
    """Per-pool targets of the pools of the external labels."""
    return expand("resolved/{pool}/{modality}.done", pool=split_pools(), modality=MODALITY)


def produce_targets_dynamic(wildcards) -> list:
    if "identify_swaps" in STAGES:
        return _expand_pools_donors(targets_identify_swaps()[0], wildcards)
    if "split_bams" in STAGES:
        return _expand_pools_donors(targets_SplitBams()[0], wildcards)
    raise ValueError(
        f"produce_targets_dynamic: until '{UNTIL}' has no per-donor targets"
    )
