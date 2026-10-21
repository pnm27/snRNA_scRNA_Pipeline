# TODO

Items are ordered roughly by priority: bugs first, then the config/schema
groundwork that the larger refactors depend on.

---

## Bugs to fix first

- [x] `demux_type: null` crashes: `get_inputs_demux` (`lib/io.py`), `get_params_demux` (`lib/params.py`) and `get_filt_barcodes` call `config['demux_type'].lower()` on `None`. Only the string `"none"` works today, although the config comment says `null` is allowed.
- [x] `get_inp_splitBam` (`lib/io.py`) checks `split_by` instead of `demux_type` in its hashsolo branch, so splitting by raw calico_solo output returns `None`.
- [x] `vcf_type` wildcard issue in `demux_samples`: `_inputs_vireo` / `_inputs_both` read `wildcards.vcf_type` when `last_step` ends in `multi_vcf`, but the rule's output path has no `{vcf_type}`. Add `{vcf_type}` to the output fold_struct for multi-VCF runs.
- [x] `get_vir_inputs`: with a 2-column `vcf_info` and `w_gt`, `donorFile` was set to `genome_1k_ref` (the 1kGP sites VCF); now `donors_per_pool_file`.
- [x] `get_cmd_str_vireo` reads `vcf_info` unconditionally and indexes `input['donorFile']`, which fails when `vcf_info` is null or no donor file is set.
- [x] Resource functions are defined twice with different values (`allocate_mem_cS`, `allocate_mem_vS`, etc. in both `resources.smk` and `genotype_demux.smk`; the later include wins). Also: typo `allocate_timem_RN`, and `allocate_mem_QM` takes no `attempt`. Resolved by the resources item below.
- [x] Multiome branch of `_inputs_vireo` lists the same `donors_classification` path twice (now one per modality).
- [x] `demux_type: both` targets didn't match the single h5ad written by `demux_samples`.
- [x] `STARsolo_kb_solo` target was missing the `modality` expansion.
- [x] VCF-suffix memory/time scaling for `cellSNP`/`vireoSNP` never applied (`'vcf_type' in wildcards` checks values, not names).
- [x] `new_config.example.yaml` was missing `cellranger_arc_count.ref`.
- [ ] Multi-VCF targets in `produce_targets.smk` use undefined `VCF_TYPE` and `suff_h5ad` (see Compatibility → Multi-VCF).
- [ ] Multiome `demux_samples` runs once per modality and writes the same combined h5ad twice (output path carries `{modality}`).

---

## Groundwork: config & schemas

- [ ] Write schemas and re-enable `validate()` in the `Snakefile` (`schema/config.schema.json` is currently gitignored). Consider growing `lib/config_models.py` (`DemuxConfig`) into typed config.
- [ ] Restructure `new_config.yaml` to separate:
  - [ ] static resources
  - [ ] experiment definitions
  - [ ] workflow behavior (RE-ASSESS `last_step`):

    **OLD STYLE**

    ```yaml
    demux_type: solo
    demux_run_type: wo_gt
    hto_demux_type: single
    ```

    **NEW STYLE**

    ```yaml
    experiment:
      demultiplexing:
        method: solo
        genotype:
          enabled: false

    hto:
      enabled: true
      multiplexing: single
    ```

  - [ ] execution modules
  - [ ] output naming conventions
- [ ] Move pool calculation logic from `input_processing.smk` into config (inspired by `downloadFromSynapse.py`, which is not in this repo yet).
- [ ] Add a tiny test dataset (one small pool) and `snakemake -n` dry runs for each `demux_type`.

---

## Resources

- [x] Move all resource definitions out of the rule files into `rules/resources.smk`, so the rules only reference them (goal: readability of the rules). Rules use `**get_resources("<rule>")`; values are an int or `{base, step, max}` (`base + step * (attempt - 1)`), memory in MB and `time_min` in minutes. A `resources:` block in the config overrides any rule:

  ```yaml
  resources:
    vireoSNP:
      cpus_per_task: 7
      mem_mb: {base: 30000, step: 1000, max: 60000}
      time_min: 240
  ```

  - [x] Support overrides by VCF suffix for `cellSNP` / `vireoSNP` (`by_vcf_suffix`).

---

## Demultiplexing outputs & stats (one design)

Design the output/stats file layout once; the items below follow from it.

- [ ] Split `demux_samples` into 2 rules:
  - [ ] RULE A: an h5ad is created after vireo output
  - [ ] RULE B: `demux_samples` picks the proper version of vireo demultiplexing, corrects swaps and creates a final h5ad
- [ ] Split the demultiplexing info into basic stats (common file), `hashsolo_stats` and `vireo_stats`, saved per demux per config/params. Mark the final version for each pool in `swap_correction_df`.
- [ ] When adding calico_solo or vireo results, take the existing demux stats file as input and append to it.
- [ ] Fix automatic output directory selection in `demux_samples` (`_output_cfg`).
- [ ] Clear strategy for `final_count_matrix_dir` across demux scenarios:
  - Both tools run simultaneously
  - Tools run sequentially (name each run separately or record the order)
- [ ] Replace separate demux config dirs (e.g. `wo_gt`, `w_gt`) with a suffix set at the top of the config file.
- [ ] `get_inp_splitBam` with `split_by: input: h5ad` always reads the hashsolo final h5ad, even for vireo.
- [ ] Use the `Final_swap_correction` file to produce the *Final* per-donor bams.
  *BEFORE*: for each group of samples that require the same round of demultiplexing, run with **new_config.yaml** set up like:

  ```yaml
  gt_check: yes # yes or no
  select_fastqs: subsetA_samples
  split_bams_pipeline:
    subset_chr: null
  ```

- [ ] Retain vireo rerun information in both the demultiplex info file and `update_logs`:
  - [x] *swap correction* addresses proper demultiplex file per each pool.
  - [ ] *swap correction* has appropriately identified donors per each pool.

---

## Demultiplexing

- [ ] Option to reuse a previous `cellSNP` run for `vireoSNP` (e.g. `gt_demux_pipeline.reuse_cellsnp_dir` read in `get_vir_inputs`); useful when there's a vcf for a w_gt run but the SNPs haven't changed (many WGS). Needs the `get_cmd_str_vireo` fix above.
- [ ] Option to use a per-pool vcf only for pileups and not for vireo (e.g. `pool_vcf_use: [pileup, vireo]`).
- [ ] Handle pools whose `cellSNP.base.vcf.gz` is empty. vireo fails before `create_inp_splitBams` runs, so either make cellSNP/vireo a checkpoint or have vireo write a sentinel the checkpoint reads.
- [ ] Retain variant matching stats (variants matched to donor VCF out of total pileup): add a `log:` directive to `vireoSNP` and parse "variants matched to donor VCF" in `update_logs.py`, instead of searching `.snakemake/log`.
- [ ] Fix `demultiplex_helper_funcs.py` for double HTOs in `parse_file`.
- [ ] Resolve the multi-HTO setup in `kite.smk` (document what currently fails).
- [ ] Beautify `get_filt_barcodes` in `genotype_demux.smk`.

---

## Wildcards & modality (one refactor)

Prerequisite for the Compatibility section.

- [ ] `wildcards.modality` handles cDNA/RNA, HTO or ATAC rather than functions (input, params, etc.).
- [ ] Simplify the wildcard structure so the folder structure encodes wildcards naturally.
- [ ] Handle varying output dirs for:
  - Demux output (single method vs. simultaneous both)
  - BAM splitting (finalizing vs. genotype purposes)
- [ ] Use the `branch` function to handle rules with varied inputs.
- [ ] Simplify target functions in `produce_targets` (`solo` / `both` targets are still TODO).

---

## Aligner abstraction

- [ ] Remove dependency on STARsolo as the only aligner.
- [ ] Add `cellranger count` (GEX-only); only `cellranger-arc count` exists today.
- [ ] Add input function for rules using `genefull_matrices` to select either `Gene` or `GeneFull` per project. Make the `gene_*` / `genefull_*` config paths consistent (`.mtx` vs `.mtx.gz`).
- [ ] STARsolo WASP mode: `--waspOutputMode SAMtag --varVCFfile` can already go through `extra_params`; the work is supplying a per-pool VCF and adding `vW` to `SAM_attr`. `wo_gt` runs have no VCF until after demultiplexing.

---

## Compatibility

- [ ] Split the design for sc/snRNA-seq, scATAC and multiome; the pipeline will be restructured.
- [ ] Make the pipeline compatible with:
  - [ ] Multiome
  - [ ] Multi-module
  - [ ] Multi-VCF
- [ ] Generalize outputs to accommodate multiome:
  - [ ] `split_bams.smk` (entire split bam pipeline)
  - [ ] `filt_chr_bams` (→ `filt_chr_bams_multiome`)
- [ ] Permanent fix for rule `filt_chr_bams_multiome`.
- [ ] Add `CollectInsertSizeMetrics` for the ATAC modality in multiome (`picard_metrics`).
- [ ] Add multiome support in `update_logs.py`.
- [ ] Support ArchR and Signac QC for ATAC QC metrics.

---

## run_update_logs.sh (one change)

- [x] Build `cmd_args` as a bash array, not a concatenated string (paths with spaces currently break).
- [x] Add a lightweight `check_required_args` for REQUIRED ARGS that fails fast with a clear message.
- [x] Make file suffixes easy to modify (`opt_values`: empty = default suffix, remove a line to skip that output).

- [x] `append_optional_flags` warns (not fails) if a key is not in `opt_flags`.

---

## Code Quality

- [ ] Python logging with the `logging` module in the remaining helper scripts (`update_logs.py` already uses it).
- [x] Remove `run` directive from `kite`, `pheno_demux3`, and `STARsolo`
- [x] Combine sub-workflows `split_bams` and `split_bams_gt`
- [x] Add config options for extra software parameters

---

## Experimental / Future

- [ ] Add new Picard metrics
- [ ] Incorporate git submodules for external repos (e.g. kite `featuremap`)

---

## Documentation

- [ ] Add tutorials:
  - Pooled snRNA-seq: single wildcard, multiple wildcards
  - scRNA-seq: single wildcard, multiple wildcards
  - Double HTOs
- [ ] Search ranking config for ReadTheDocs (`search.ranking` in `.readthedocs.yaml`)

---

## Benched

- [ ] Dash app: pair-wise sample comparison (donor vs donor or pool vs pool) sub-plot pop up.
- [ ] Agentic AI support for genotype-based demultiplexing using Ollama (with support for public vLLMs).
- [ ] Shannon entropy for a better selection of SNPs (benchmark against the existing `max50k` / `25pct` sets first).
