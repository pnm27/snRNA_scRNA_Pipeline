#!/usr/bin/env python3


from dataclasses import dataclass, field
from typing import List


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
        demux_type = config.get("demux_type")
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
                    fold_struct=config['fold_struct_gt_demux'].rstrip('/')
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
                    fold_struct=config['fold_struct_gt_demux'].rstrip('/')

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

        else:
            raise ValueError(
                f"Unsupported demux_type: {demux_type}"
            )

        return cls(
            demuxes=demuxes,
            multiome="multiome" in config["last_step"].lower(),
            create_h5ad=create_h5ad
        )