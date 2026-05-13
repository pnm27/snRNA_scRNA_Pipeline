from lib.io import (
    get_filt_barcodes, 
    get_cellsnp_inputs,
    get_vir_inputs,
    get_inputs_demux
)
from lib.params import get_params_demux
from functools import partial


# ── Output resolver ───────────────────────────────────────────────────────────

def _output_cfg(config):
    """Return (matrix_dir, fold_struct, matrix_h5ad, demux_info_dir, demux_info, suffix)."""
    demux_type = config['demux_type'].lower()

    if demux_type in ('solo', 'add_solo'):
        c = config['hashsolo_demux_pipeline']
        return (c['final_count_matrix_dir'], config['fold_struct_demux'],
                c['final_count_matrix_h5ad'],
                c['demultiplex_info_dir'], c['demultiplex_info'], "")

    if demux_type in ('vireo', 'add_vireo'):
        c = config['gt_demux_pipeline']
        return (c['final_count_matrix_dir'], config['fold_struct_demux'],
                c['final_count_matrix_h5ad'],
                c['demultiplex_info_dir'], c['demultiplex_info'], "")

    if demux_type == 'both':
        c = config['demultiplex']
        return (c['demux_count_matrix_dir'], config['fold_struct_demux'],
                c['final_count_matrix_h5ad'],
                c['demux_count_matrix_dir'], c['demultiplex_info'], "")

    # create_h5ad_only — strip the "_vS" suffix
    c = config['gt_demux_pipeline']
    return (c['final_count_matrix_dir'], config['fold_struct'],
            c['final_count_matrix_h5ad'].replace("_vS", ""),
            c['demultiplex_info_dir'],
            c['demultiplex_info'].replace("_vS", ""), "")


def _make_outputs(config):
    mat_dir, fold, mat_h5ad, demux_dir, demux_info, _ = _output_cfg(config)
    return (
        f"{mat_dir}{fold}{mat_h5ad}",
        f"{demux_dir}{fold}{demux_info}",
    )


# ── Consolidated rule ─────────────────────────────────────────────────────────

_out = _make_outputs(config)

rule demux_samples:
    input:
        lambda wildcards: get_inputs_demux(wildcards, config)

    output:
        _out[0],
        _out[1]

    params:
        pool_name = lambda wildcards: wildcards.pool,
        extra     = partial(get_params_demux, config=config)

    resources:
        mem_mb   = lambda wildcards, attempt: 3500 * attempt + 3500,
        time_min = lambda wildcards, attempt: 15   * attempt + 15

    conda: "../envs/basic_sctools.yaml"

    shell:
        """
        trap 'echo "Received SIGTERM, killing children"; kill 0; wait' SIGTERM
        python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
        """
