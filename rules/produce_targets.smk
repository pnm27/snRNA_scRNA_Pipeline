import pandas as pd
import os
from typing import Literal


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# STATIC, pool-only expansion
def _expand_pools(template: str) -> list:
    """Expand a path template over all pools."""
    return expand(template, pool=POOLS)


# DYNAMIC, pool+donors expansion
def _expand_pools_donors(template: str, wildcards) -> list:
    """Expand a path template over pool + its donors, zipped."""
    wc_dict = get_real_donors({"pool": wildcards.pool})
    return expand(template, zip, **wc_dict)


def _picard_targets(progs: Literal["all", "gc", "rnaseq"] = "all") -> list:
    """Return Picard output paths for the requested program(s)."""
    inp_pref = config["STARsolo_pipeline"]["bams_dir"]
    f_st = config["fold_struct"]
    files = {
        "rnaseq": config["picard_pipeline"]["rnaseq_metrics"],
        "gc":     config["picard_pipeline"]["gc_summary_metrics"],
    }
    keys = files.keys() if progs == "all" else [progs]
    return [f"{inp_pref}{f_st}{files[k]}" for k in keys]


# ---------------------------------------------------------------------------
# Per-step target builders
# ---------------------------------------------------------------------------

def targets_STARsolo() -> list:
    inp_pref = config["STARsolo_pipeline"]["bams_dir"]
    bai_suff = config["STARsolo_pipeline"]["bai"]
    f_st     = config["fold_struct"]
    return [f"{inp_pref}{f_st}{bai_suff}"]


def targets_PICARD(progs: Literal["all", "gc", "rnaseq"] = "all") -> list:
    return _picard_targets(progs)


def targets_cellSNP(multiome: bool = False) -> list:
    sub_dirs  = ["ATAC", "cDNA"] if multiome else [""]
    out_dir   = config["gt_demux_pipeline"]["cellsnp_dir"]
    fs_gt     = config["fold_struct_gt_demux"]
    cells_suff = config["gt_demux_pipeline"]["cellsnp_cells"]
    return [
        os.path.join(f"{out_dir}{fs_gt}", d, cells_suff)
        for d in sub_dirs
    ]


def targets_demux(suffix: list, out_dir: list | None,
    h5ad: bool = False,
    # demux_type: Literal['vireo', 'solo', 'both'],
    multiome: bool = False) -> list:
    sub_dirs = ["ATAC", "cDNA"] if multiome else [""]
    ret_targets = []
    if not h5ad:
        # NEED TO DO FOR SOLO
        for suff in suffix:
            out_dir = config["gt_demux_pipeline"]["vireosnp_dir"]
            fs_gt   = config["fold_struct_gt_demux"]
            suff    = config["gt_demux_pipeline"]["donors_classification"]
            ret_targets.append(
                os.path.join(f"{out_dir}{fs_gt}", d, suff)
                for d in sub_dirs
            )
    else:
        for od, suff in zip(out_dir, suffix):
            #  NOT SAME FOR CREATE_H5AD_ONLY
            fs_gt   = config["fold_struct_demux"]
            # out_dir = config["gt_demux_pipeline"]["final_count_matrix_dir"]
            # suff    = config["gt_demux_pipeline"]["final_count_matrix_h5ad"]
            ret_targets.append(os.path.join(f"{od}{fs_gt}{suff}"))

    return ret_targets


def targets_gt_demux2(multiome: bool = False) -> list:
    sub_dirs = ["ATAC", "cDNA"] if multiome else [""]
    out_dir  = config["gt_demux_pipeline"]["final_count_matrix_dir"]
    fs_gt    = config["fold_struct_demux"]
    return [
        os.path.join(out_dir, d, fs_gt)
        for d in sub_dirs
    ]


def targets_SplitBams(multiome: bool = False) -> list:
    sub_dirs = ["ATAC", "cDNA"] if multiome else [""]
    out_dir  = (
        config["split_bams_pipeline"]["split_bams_proxy_dir2"]
        if config["gt_check"]
        else config["split_bams_pipeline"]["split_bams_proxy_dir"]
    )
    op = config["fold_struct_bam_split1"]
    return [
        os.path.join(out_dir, op, d) if d else os.path.join(out_dir, op)
        for d in sub_dirs
    ]


def targets_gt_demux_identify_swaps(multiome: bool = False) -> list:
    out_dir = config["identify_swaps"]["mbv_out_dir"]
    fs      = config["fold_struct_swaps_check"]
    suff    = config["identify_swaps"]["mbv_suffix"]
    return [os.path.join(out_dir, f"{fs}{suff}")]


def targets_multiome(last: str, suffix: list | None, 
    out_dir: list | None,
    h5ad: bool = False,
    multiome: bool = False) -> list:
    # required_args = ['vireo', 'solo']
    dispatch = {
        "alignment":   lambda: [
            f"{config['cellranger_arc_count']['bams_dir']}{{pool}}/{f}"
            for f in [
                "filtered_feature_bc_matrix/barcodes.tsv.gz",
                "filtered_feature_bc_matrix/features.tsv.gz",
                "filtered_feature_bc_matrix/matrix.mtx.gz",
            ]
        ],
        "vireo": lambda x, y: targets_demux(x, y, h5ad=h5ad, multiome=multiome),
        "splitBams":   lambda: targets_SplitBams(multiome=multiome),
        "identifySwaps": lambda: targets_gt_demux_identify_swaps(multiome=multiome),
        "solo": [], # ADD LATER
        "both": [], # ADD LATER
    }
    if last not in dispatch:
        raise ValueError(f"Unknown multiome target stage: '{last}'")

    return dispatch[last](suffix, out_dir) if suffix is not None else dispatch[last]()


# ---------------------------------------------------------------------------
# produce_targets  — static expansion (DAG construction time)
# ---------------------------------------------------------------------------

def produce_targets(wildcards) -> list:
    step     = config["last_step"].lower()
    metrics  = (config["picard_metrics"].lower()
                if config["picard_metrics"] is not None else None)
    h5ad     = config.get("create_final_h5ad", False) if "gt_demux" in step else False
    demux_type = config.get("demux_type", None)
    h5ad_suffixes = []
    h5ad_outdir = []
    multiome = "multiome" in step
    suff_h5ad = config["gt_demux_pipeline"]["final_count_matrix_h5ad"]


    if not multiome and demux_type in ['vireo', 'add_vireo']:
        demux_type = 'vireo'
        h5ad_outdir = [
            config['gt_demux_pipeline']['final_count_matrix_dir'],
            ]
        h5ad_suffixes = [
            config['gt_demux_pipeline']['final_count_matrix_h5ad'],
            ]
    elif not multiome and demux_type in ['solo', 'add_solo']:
        demux_type = 'solo'
        h5ad_outdir = [
            config['hashsolo_demux_pipeline']['final_count_matrix_dir'],
            ]
        h5ad_suffixes = [
            config['hashsolo_demux_pipeline']['final_count_matrix_h5ad'],
            ]
    #  TO DO
    elif not multiome and demux_type == 'both':
        h5ad_suffixes = [
            config['gt_demux_pipeline']['final_count_matrix_h5ad'],
            config['hashsolo_demux_pipeline']['final_count_matrix_h5ad'],
            ]

    def _with_picard(targets: list) -> list:
        """Optionally append Picard targets to an existing list."""
        if metrics is not None:
            targets += _expand_pools(_picard_targets(progs=metrics))
        return targets

    # -- dispatch -----------------------------------------------------------
    if step == "all_multi_vcf":
        base = [t for t in _expand_pools(targets_STARsolo()[0])]
        gt2  = targets_gt_demux2()
        if isinstance(VCF_TYPE, list):
            vcf_targets = [
                expand(f"{t}_{{vcf_type}}{suff_h5ad}", vcf_type=VCF_TYPE)
                for t in _expand_pools(gt2[0])
            ]
            return base + [f for sub in vcf_targets for f in sub]
        elif isinstance(VCF_TYPE, str):
            return base + [
                expand(f"{t}_{VCF_TYPE}{suff_h5ad}", pool=POOLS)
                for t in gt2
            ]
        else:
            raise ValueError("VCF_TYPE must be a str or list.")

    elif step == "starsolo_kb_solo":
        return _expand_pools(targets_STARsolo()[0])

    elif step == "starsolo_cellsnp":
        return _with_picard(_expand_pools(targets_cellSNP(multiome=False)[0]))

    elif step == "starsolo_demux":
        fs_gt = config['gt_demux_pipeline']['final_count_matrix_h5ad']
        return _with_picard(
            _expand_pools(targets_demux(fs_gt, h5ad=h5ad,
            multiome=multiome)[0])
        )

    elif step == "starsolo_gt_demux_multi_vcf":
        gt2 = targets_gt_demux2()
        if isinstance(VCF_TYPE, list):
            return [
                f
                for t in _expand_pools(gt2[0])
                for f in expand(f"{t}_{{vcf_type}}{suff_h5ad}", vcf_type=VCF_TYPE)
            ]
        elif isinstance(VCF_TYPE, str):
            return expand(f"{gt2[0]}_{VCF_TYPE}{suff_h5ad}", pool=POOLS)
        else:
            raise ValueError("VCF_TYPE must be a str or list.")

    elif step == "multiome_alignment":
        return [f for t in targets_multiome("alignment", None)
                for f in _expand_pools(t)]

    elif step == "multiome_gt_demux":
        h5ad_outdir = [config['gt_demux_pipeline']['final_count_matrix_dir']]
        h5ad_suffixes = [config['gt_demux_pipeline']['final_count_matrix_h5ad']]
        return [
            f for t in targets_multiome("vireo", h5ad_suffixes, h5ad_outdir, h5ad=h5ad, multiome=multiome)
            for f in _expand_pools(t)
        ]

    elif step == "multiome_split_bams_gt_demux":
        targets = targets_multiome("splitBams", None, h5ad=h5ad, multiome=multiome)
        return [
            f
            for i, t in enumerate(targets)
            for f in _expand_pools(f"{t}.txt" if i <= 1 else t)
        ]

    elif step == "multiome_gt_demux_identify_swaps":
        return [f for t in targets_multiome("identifySwaps", None, multiome=multiome)
                for f in _expand_pools(t)]

    elif "starsolo" in step:
        ret_list=[]
        if metrics is None:
            ret_list += _expand_pools(targets_STARsolo()[0])
        else:
            ret_list += _expand_pools(_picard_targets(progs=metrics))

        if any(s in step for s in ["split_bams", "identify_swaps"]):
            ret_list += expand("resolved/{pool}.done", pool=POOLS)

        if h5ad:
            ret_list += _expand_pools(
                targets_demux(h5ad_suffixes, h5ad_outdir, h5ad=h5ad, multiome=multiome)[0])
        return ret_list


    else:
        raise ValueError(f"Unrecognised last_step: '{step}'")


# ---------------------------------------------------------------------------
# produce_targets_dynamic — runtime expansion (after hash file exists)
# ---------------------------------------------------------------------------

def produce_targets_dynamic(wildcards) -> list:
    step = config["last_step"].lower()

    if step == "starsolo_split_bams":
        return _expand_pools_donors(targets_SplitBams()[0], wildcards)

    elif step == "starsolo_split_bams_gt_demux":
        targets = targets_SplitBams()
        return [
            f
            for i, t in enumerate(targets)
            for f in _expand_pools_donors(f"{t}.txt" if i == 0 else t, wildcards)
        ]

    elif step == "starsolo_split_bams_gt_demux_multi_vcf":
        suff_h5ad = config["gt_demux_pipeline"]["final_count_matrix_h5ad"]
        targets   = targets_SplitBams()
        if isinstance(VCF_TYPE, list):
            return [
                f
                for t in _expand_pools_donors(targets[0], wildcards)
                for f in expand(f"{t}_{{vcf_type}}{suff_h5ad}", vcf_type=VCF_TYPE)
            ]
        elif isinstance(VCF_TYPE, str):
            return [
                f"{t}_{VCF_TYPE}{suff_h5ad}"
                for t in _expand_pools_donors(targets[0], wildcards)
            ]
        else:
            raise ValueError("VCF_TYPE must be a str or list.")

    elif step == "starsolo_gt_demux_identify_swaps":
        return _expand_pools_donors(
            targets_gt_demux_identify_swaps(multiome=False)[0], wildcards
        )

    else:
        raise ValueError(
            f"produce_targets_dynamic: unrecognised last_step '{step}'"
        )
        