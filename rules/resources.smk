# Resources for every rule, kept here so the rules stay readable.
# Memory is in MB and time is in minutes.
#
# Each value can be:
#   - an int, used as is
#   - a dict, scaled with the retry attempt and capped at 'max' (optional):
#       base + step * (attempt - 1)
#   - a function (wildcards, attempt), for special cases
#
# 'by_vcf_suffix' (optional) overrides a dict value when wildcards.vcf_type
# ends with one of its keys (checked in order). The key "" matches any
# vcf_type, so it is the fallback for runs that have one.
#
# Values under 'resources' in the config (same layout) override these, e.g.
#   resources:
#     vireoSNP:
#       mem_mb: {base: 30000, step: 1000}
#
# Usage in a rule:
#   resources:
#       **get_resources("cellSNP")


# Per-pool max RAM for STAR comes from the optional second column of
# 'select_fastqs' (SS_MEM, see input_processing.smk)
def allocate_mem_SS(wildcards, attempt):
    if not SS_MEM or wildcards.pool not in SS_MEM:
        return 45000+1000*(attempt-1) # DEFAULT
    else:
        return SS_MEM[wildcards.pool] +4000*(attempt-1)


RESOURCES = {
    # STARsolo.smk
    "STARsolo_sort": {
        "cpus_per_task": 10,
        "mem_mb": allocate_mem_SS,
        "time_min": 1440,
    },
    "samtools_index": {
        "cpus_per_task": 4,
        "mem_mb": 5000,
        "time_min": 720,
    },

    # picard_metrics.smk
    "Picard_GC_bias_metrics": {
        "cpus_per_task": 2,
        "mem_mb": {"base": 7500, "step": 2500},
        "time_min": {"base": 360, "step": 60},
    },
    "Picard_RNAseq_metrics": {
        "cpus_per_task": 4,
        "mem_mb": {"base": 2500, "step": 2500},
        "time_min": {"base": 210, "step": 30},
    },

    # kite.smk
    "create_FB": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 130, "step": 50},
        "time_min": {"base": 10, "step": 5},
    },
    "create_mismatch_fasta": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 120, "step": 50},
        "time_min": {"base": 10, "step": 5},
    },
    "build_kallisto_index": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 2000, "step": 1000},
        "time_min": {"base": 10, "step": 5},
    },
    "run_kallisto": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 1500, "step": 1000},
        "time_min": {"base": 40, "step": 20},
    },
    "run_bustools_correct": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 3000, "step": 1500},
        "time_min": {"base": 10, "step": 5},
    },
    "run_bustools_sort": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 6000, "step": 3000},
        "time_min": {"base": 10, "step": 5},
    },
    "run_bustools_count": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 3000, "step": 1500},
        "time_min": {"base": 10, "step": 5},
    },
    "create_h5ad_bustools": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 3000, "step": 1500},
        "time_min": {"base": 10, "step": 5},
    },

    # calico_solo_demux.smk
    "run_calico_solo": {
        "cpus_per_task": 2,
        "mem_mb": {"base": 7000, "step": 4500},
        "time_min": {"base": 10, "step": 5},
    },

    # genotype_demux.smk
    "create_inp_cellSNP": {
        "cpus_per_task": 2,
        "mem_mb": {"base": 2000, "step": 500},
        "time_min": {"base": 3, "step": 2},
    },
    "cellSNP": {
        "cpus_per_task": 8,
        "mem_mb": {
            "base": 200, "step": 200,
            "by_vcf_suffix": {
                "_max50k": {"base": 200, "step": 80},
                "25pct": {"base": 800, "step": 100},
                "gencode": {"base": 1300, "step": 200},
                "": {"base": 1500, "step": 200},
            },
        },
        "time_min": {
            "base": 240, "step": 60,
            "by_vcf_suffix": {
                "_max50k": {"base": 30, "step": 20},
                "25pct": {"base": 70, "step": 20},
                "gencode": {"base": 180, "step": 40},
            },
        },
    },
    "vireoSNP": {
        "cpus_per_task": 7,
        "mem_mb": {
            "base": 35000, "step": 500,
            "by_vcf_suffix": {
                "_max50k": {"base": 3000, "step": 100},
                "25pct": {"base": 5000, "step": 400},
                "gencode": {"base": 10000, "step": 600},
                "": {"base": 22000, "step": 500},
            },
        },
        "time_min": {
            "base": 180, "step": 60,
            "by_vcf_suffix": {
                "_max50k": {"base": 10, "step": 5},
                "25pct": {"base": 70, "step": 20},
                "gencode": {"base": 180, "step": 40},
            },
        },
    },

    # demultiplex.smk
    "demux_run": {
        "mem_mb": {"base": 7000, "step": 3500},
        "time_min": {"base": 30, "step": 15},
    },
    "demux_final": {
        "mem_mb": {"base": 7000, "step": 3500},
        "time_min": {"base": 15, "step": 15},
    },

    # split_bams.smk
    "create_inp_splitBams": {
        "cpus_per_task": 2,
        "mem_mb": 1000,
        "time_min": 10,
    },
    "bamfilt_by_CB": {
        "cpus_per_task": 1,
        "mem_mb": 250,
        "time_min": 300,
    },
    "filt_chr_bams": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 1000, "step": 200},
        "time_min": {"base": 20, "step": 20},
    },
    "filt_chr_bams_multiome_cdna": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 1000, "step": 200},
        "time_min": {"base": 20, "step": 20},
    },
    "filt_chr_bams_multiome_atac": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 1000, "step": 200},
        "time_min": {"base": 20, "step": 20},
    },
    "split_bams": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 250, "step": 100},
        "time_min": 30,
    },

    # identify_swaps.smk
    "qtltools_mbv": {
        "cpus_per_task": 1,
        "mem_mb": 5000,
        "time_min": {"base": 60, "step": 20},
    },

    # deeptools.smk
    "multiBamSummary": {
        "cpus_per_task": 5,
        "mem_mb": {"base": 400, "step": 200},
        "time_min": {"base": 300, "step": 30},
    },
    "plotCorrelation": {
        "cpus_per_task": 2,
        "mem_mb": {"base": 200, "step": 500},
        "time_min": {"base": 10, "step": 200},
    },

    # cellranger.smk
    "inp_cellranger_arc_count": {
        "cpus_per_task": 1,
        "mem_mb": {"base": 2000, "step": 1000},
        "time_min": 10,
    },
    "cellranger_arc_count": {
        "cpus_per_task": 20,
        "mem_mb": {"base": 16000, "step": 2000},
        "time_min": 2880,
    },
}


def _resolve(spec, wildcards, attempt):
    """Value of one resource spec for a job (see the top of this file)."""
    if callable(spec):
        return spec(wildcards, attempt)

    vcf_type = getattr(wildcards, 'vcf_type', None)
    if vcf_type is not None:
        for suffix, sub_spec in spec.get("by_vcf_suffix", {}).items():
            if vcf_type.endswith(suffix):
                spec = sub_spec
                break

    value = spec["base"] + spec.get("step", 0) * (attempt - 1)
    return min(value, spec.get("max", value))


def get_resources(rule):
    """Keyword arguments for the 'resources' directive of a rule."""
    specs = {
        **RESOURCES[rule],
        **((config.get("resources") or {}).get(rule) or {}),
    }

    ret = {}
    for name, spec in specs.items():
        if isinstance(spec, int):
            ret[name] = spec
        else:
            ret[name] = (
                lambda wildcards, attempt, spec=spec:
                    _resolve(spec, wildcards, attempt)
            )
    return ret
