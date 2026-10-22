from lib.io import (
    get_inputs_demux_run,
    get_inputs_demux_final,
)
from lib.params import get_params_demux, get_params_demux_final
from lib.config_models import demux_tools, demux_run, final_paths
from functools import partial


# Output paths of both rules are in lib/config_models.py

# ── Rule A: one h5ad per demultiplexing tool of this run ─────────────────────

for _tool in demux_tools(config):
    _run = demux_run(config, _tool, config['experiment']['demultiplexing']['run_tag'])

    rule:
        name: f"demux_run_{_tool}"

        input:
            unpack(partial(get_inputs_demux_run, config=config, tool=_tool))

        output:
            h5ad = _run.h5ad,
            **({'stats': _run.stats} if _run.stats else {})

        params:
            pool_name = lambda wildcards: wildcards.pool,
            extra     = partial(get_params_demux, config=config, tool=_tool)

        resources:
            **get_resources("demux_run")

        conda: "../envs/basic_sctools.yaml"

        shell:
            """
            trap 'echo "Received SIGTERM, killing children"; kill 0; wait' SIGTERM
            python3 helper_py_scripts/demul_samples.py {params.extra} --pool_name {params.pool_name}
            """


# ── Rule B: final h5ad of a pool ─────────────────────────────────────────────
# Combines the Rule A outputs picked for the pool: its row in the swap
# correction file, else the tools of this run (plus, for add_solo/add_vireo,
# the other tool from 'previous_run_tag')

rule demux_final:
    input:
        unpack(partial(get_inputs_demux_final, config=config))

    output:
        h5ad        = final_paths(config).h5ad,
        basic_stats = final_paths(config).basic_stats,
        info        = final_paths(config).info

    params:
        extra = partial(get_params_demux_final, config=config)

    resources:
        **get_resources("demux_final")

    conda: "../envs/basic_sctools.yaml"

    shell:
        """
        python3 helper_py_scripts/demux_final.py --h5ad {input.h5ad} \
            --stats {input.stats} --pool_name {wildcards.pool} \
            -o {output.h5ad} --basic_stats {output.basic_stats} \
            --info_out {output.info} {params.extra}
        """


# ── Demultiplex versions used for the final h5ad of each pool ────────────────

rule collect_demux_versions:
    input:
        # Multi-VCF runs: one table per vcf_type
        expand(final_paths(config).basic_stats, pool=POOLS, allow_missing=True)

    output:
        final_paths(config).versions

    params:
        pools = POOLS

    resources:
        **get_resources("collect_demux_versions")

    shell:
        """
        python3 helper_py_scripts/collect_demux_versions.py \
            --pools {params.pools} --basic_stats {input} -o {output}
        """
