#!/usr/bin/env python3


from dataclasses import dataclass, field
from typing import List
from .utils import get_demux_type


@dataclass
class DemuxResult:
    demux_type: str
    outdir: str
    suffix: str
    fold_struct: str


@dataclass
class DemuxConfig:

    demuxes: list[DemuxResult] = field(default_factory=list)
    multiome: bool = False
    create_h5ad: bool = False

    @classmethod
    def from_snakemake(cls, config):

        demuxes = []
        demux_type = get_demux_type(config)
        create_h5ad = config.get(
                "create_final_h5ad",
                False
            )
        if demux_type in ["vireo", "add_vireo"]:
            demuxes.append(
                DemuxResult(
                    demux_type="vireo",
                    outdir=(
                        config['gt_demux_pipeline']['final_count_matrix_dir']
                        if create_h5ad
                        else config["gt_demux_pipeline"]["vireosnp_dir"]
                    ),
                    suffix=(
                        config['gt_demux_pipeline']['final_count_matrix_h5ad']
                        if create_h5ad
                        else config["gt_demux_pipeline"]["donors_classification"]
                    ),
                    fold_struct=(
                        config['fold_struct_gt_demux'].rstrip('/')
                        if create_h5ad
                        else config['fold_struct_gt_demux']
                    )
                )
            )

        elif demux_type in ["solo", "add_solo"]:
            demuxes.append(
                DemuxResult(
                    demux_type="solo",
                    outdir=config[
                        'hashsolo_demux_pipeline'
                    ]['final_count_matrix_dir'],
                    suffix=config[
                        'hashsolo_demux_pipeline'
                    ]['final_count_matrix_h5ad'],
                    fold_struct=config['fold_struct_demux']
                )
            )

        elif demux_type == "both" and create_h5ad:
            # Both demuxes are combined into one h5ad by 'demux_samples'
            # (see '_output_cfg' in rules/demultiplex.smk)
            demuxes.append(
                DemuxResult(
                    demux_type="both",
                    outdir=config[
                        'demultiplex'
                    ]['demux_count_matrix_dir'],
                    suffix=config[
                        'demultiplex'
                    ]['final_count_matrix_h5ad'],
                    fold_struct=config['fold_struct_demux']
                )
            )

        elif demux_type == "both":
            demuxes.extend([
                DemuxResult(
                    demux_type="vireo",
                    outdir=(
                        config['gt_demux_pipeline']['final_count_matrix_dir']
                        if create_h5ad
                        else config["gt_demux_pipeline"]["vireosnp_dir"]
                    ),
                    suffix=(
                        config['gt_demux_pipeline']['final_count_matrix_h5ad']
                        if create_h5ad
                        else config["gt_demux_pipeline"]["donors_classification"]
                    ),
                    fold_struct=(
                        config['fold_struct_gt_demux'].rstrip('/')
                        if create_h5ad
                        else config['fold_struct_gt_demux']
                    )
                ),
                DemuxResult(
                    demux_type="solo",
                    outdir=config[
                        'hashsolo_demux_pipeline'
                    ]['final_count_matrix_dir'],
                    suffix=config[
                        'hashsolo_demux_pipeline'
                    ]['final_count_matrix_h5ad'],
                    fold_struct=config['fold_struct_demux']
                )
            ])
        elif demux_type is None:
            # create_h5ad_only — must match the 'demux_samples' output
            # (see '_output_cfg' in rules/demultiplex.smk)
            demuxes.append(
                DemuxResult(
                    demux_type="none",
                    outdir=config[
                        'gt_demux_pipeline'
                    ]['final_count_matrix_dir'],
                    suffix=config[
                        'gt_demux_pipeline'
                    ]['final_count_matrix_h5ad'].replace("_vS", ""),
                    fold_struct=config['fold_struct']
                )
            )
        else:
            raise ValueError(
                f"Unsupported demux_type: {demux_type}"
            )

        return cls(
            demuxes=demuxes,
            multiome="multiome" in config["last_step"].lower(),
            create_h5ad=create_h5ad
        )