import os
from lib.config_models import DemuxConfig, vcf_types
from lib.aligner import bam_paths, get_aligner, out_prefix
from lib.aligner import count_matrix as aligner_count_matrix


# Targets of 'rule all' come from the stages that 'last_step' runs:
#   align           alignment (or Picard's metrics, when 'picard_metrics' is set)
#   picard          Picard's metrics, when 'picard_metrics' is set
#   cellsnp         cellSNP's pileup
#   demux           demultiplexing (see DemuxConfig)
#   split_bams      per-donor bams      } known once the donors of each pool
#   identify_swaps  QTLtools-mbv        } are (checkpoint 'create_inp_splitBams')

# last_step values whose stages don't follow from their name (see last_step_stages)
SPECIAL_STEPS = {
    "starsolo_kb_solo":            ["align"],
    "starsolo_cellsnp":            ["cellsnp", "picard"],
    "starsolo_demux":              ["demux", "picard"],
    "starsolo_gt_demux_multi_vcf": ["demux"],
    "all_multi_vcf":               ["align", "demux"],
    "cellranger_cellsnp":          ["cellsnp"],
    "cellranger_demux":            ["demux"],
    "multiome_alignment":          ["align"],
    "multiome_gt_demux":           ["demux"],
}


def last_step_stages(step: str) -> list:
    """Stages run for 'last_step' (lowercase)."""
    if step in SPECIAL_STEPS:
        return SPECIAL_STEPS[step]
    if not step.startswith(("starsolo", "cellranger", "multiome")):
        raise ValueError(f"Unrecognised last_step: '{step}'")

    stages = ["align"]
    if any(s in step for s in ["gt_demux", "split_bams", "identify_swaps"]):
        stages.append("demux")
    if any(s in step for s in ["split_bams", "identify_swaps"]):
        stages.append("split_bams")
    if "identify_swaps" in step:
        stages.append("identify_swaps")
    return stages


# ---------------------------------------------------------------------------
# This run
# ---------------------------------------------------------------------------
STEP = config["last_step"].lower()
STAGES = last_step_stages(STEP)
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
    progs = (config.get("picard_metrics") or "").lower()
    if not progs or get_aligner(config) != "starsolo":
        return []
    files = {
        "rnaseq": config["picard_pipeline"]["rnaseq_metrics"],
        "gc":     config["picard_pipeline"]["gc_summary_metrics"],
    }
    keys = files.keys() if progs == "all" else [progs]
    return [f"{out_prefix(config)}{files[k]}" for k in keys]


# ---------------------------------------------------------------------------
# Per-stage target builders (path templates)
# ---------------------------------------------------------------------------

def targets_align() -> list:
    if get_aligner(config) == "cellranger_arc":
        out_dir = config['cellranger_arc_count']['bams_dir']
        return [
            f"{out_dir}{{pool}}/filtered_feature_bc_matrix/{f}"
            for f in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]
        ]
    if get_aligner(config) == "cellranger":
        return [aligner_count_matrix(config), bam_paths(config, "cDNA")["bai"]]
    # Picard's outputs stand for the alignment when they are asked for
    return _picard_targets() or [bam_paths(config, "cDNA")["bai"]]


def targets_cellSNP() -> list:
    c = config["gt_demux_pipeline"]
    return [f"{c['cellsnp_dir']}{config['fold_struct_gt_demux']}{c['cellsnp_cells']}"]


def targets_demux() -> list:
    """Demultiplexing targets of this run (see DemuxConfig): the final h5ad
    of each pool (and the table of demultiplex versions used) or, without
    'create_final_h5ad', the demultiplexing tools' own outputs."""
    return list(CFG.targets) + list(CFG.summary_targets)


def targets_SplitBams() -> list:
    """Per-donor bams (split_bams_dir has the gt_check/ or final/ sub-folder)."""
    out_dir  = f"{config['split_bams_pipeline']['split_bams_dir']}"
    op = config["fold_struct_bam_split2"]
    return [
        os.path.join(out_dir, op, "{donor}.bam")
    ]


def targets_identify_swaps() -> list:
    out_dir = config["identify_swaps"]["mbv_out_dir"]
    fs      = config["fold_struct_swaps_check"]
    suff    = config["identify_swaps"]["mbv_suffix"]
    return [os.path.join(out_dir, f"{fs}{suff}")]


STATIC_TARGETS = {
    "align":   targets_align,
    "picard":  _picard_targets,
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
    if "split_bams" in STAGES:
        targets += expand("resolved/{pool}/{modality}.done", pool=POOLS, modality=MODALITY)

    return targets


# ---------------------------------------------------------------------------
# produce_targets_dynamic — runtime expansion (after hash file exists)
# ---------------------------------------------------------------------------

def produce_targets_dynamic(wildcards) -> list:
    if "identify_swaps" in STAGES:
        return _expand_pools_donors(targets_identify_swaps()[0], wildcards)
    if "split_bams" in STAGES:
        return _expand_pools_donors(targets_SplitBams()[0], wildcards)
    raise ValueError(
        f"produce_targets_dynamic: last_step '{STEP}' has no per-donor targets"
    )
