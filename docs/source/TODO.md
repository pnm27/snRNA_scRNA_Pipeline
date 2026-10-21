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
- [x] Multi-VCF targets in `produce_targets.smk` used undefined `VCF_TYPE` and `suff_h5ad`; `vcf_type` now comes from `vcf_info_columns: vcf`.
- [x] Multiome `demux_samples` ran once per modality and wrote the same combined h5ad twice (output path carried `{modality}`); `demux_run_vireo` / `demux_final` outputs have no `{modality}`.
- [x] `STARsolo_cellsnp` target had an extra `cDNA/` folder that no rule makes.
- [x] Split bams targets (`starsolo_split_bams`, `multiome_split_bams_gt_demux`) pointed at `split_bams_proxy_dir*` files that no rule makes.
- [x] `identify_swaps` only had a target for the first donor of each pool.
- [x] `filt_chr_bams_multiome_*` used `SimpleNamespace` without importing it.
- [x] `bamfilt_by_CB` wrote every pool's barcodes to one file literally named `{pool}_bc.txt`.

---

## Groundwork: config & schemas

- [x] Write a config schema (`schema/config.schema.yaml`) and re-enable `validate()` in the `Snakefile`; defaults fill in optional keys missing from older configs.
  - [ ] Tighten the schema once the config is restructured (required nested keys, enums for `last_step`).
  - [ ] Consider growing `lib/config_models.py` (`DemuxConfig`) into typed config.
- [x] Add a tiny test dataset (one small pool) and `snakemake -n` dry runs for each `demux_type` (`tests/dryrun/run_dryruns.sh`).
  - [x] Add scenarios for multi-VCF and the split-bams/identify-swaps steps (the harness marks the jobs up to the donor checkpoint as done, `# donors:`).

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

The layout is described in `docs/source/sub_workflows.md` (Demultiplexing outputs).

- [x] Split `demux_samples` into 2 rules:
  - [x] RULE A: `demux_run_<tool>` creates an h5ad per demultiplexing tool (and run)
  - [x] RULE B: `demux_final` picks the version per pool (swap correction file, else this run), corrects swaps (`Donor_final`) and creates the final h5ad
- [x] Split the demultiplexing info into basic stats (`_basic_stats.tsv`), `_hashsolo_stats.tsv` and `_vireo_stats.tsv`, saved per tool per run. The final version of each pool is listed in `final_demux_versions.tsv` (the swap correction file is left as it is).
- [x] calico_solo / vireo results added later (`add_solo` / `add_vireo`) are combined by `demux_final` with the other tool's run (`demultiplex.previous_run_tag`), instead of appending to its stats file.
- [x] Removed the automatic output directory selection in `demux_samples` (`_output_cfg`); paths are in `lib/config_models.py`.
- [x] Clear strategy for `final_count_matrix_dir` across demux scenarios: each tool writes to its own dir (per `run_tag`); `demultiplex.demux_count_matrix_dir` has the final h5ad of each pool.
  - Both tools run simultaneously: `demux_final` combines both runs
  - Tools run sequentially: each run has its own `run_tag`; the order is in the `Final demultiplex versions` row of the basic stats
- [x] Replaced separate demux config dirs (e.g. `wo_gt`, `w_gt`) with `run_tag` at the top of the config file.
- [x] `get_inp_splitBam` with `split_by: input: h5ad` always read the hashsolo final h5ad, even for vireo.
- [x] Use the swap correction file to produce the *Final* per-donor bams (`split_by: input: final`).
- [x] Retain vireo rerun information in both the demultiplex info file and `update_logs`:
  - [x] *swap correction* addresses proper demultiplex file per each pool.
  - [x] *swap correction* has appropriately identified donors per each pool (`Donor_final`; donor count in `final_demux_versions.tsv`).
- [ ] `update_logs.py`: read the basic/hashsolo/vireo stats files separately (its `Basic_Demultiplex_stats` etc. placeholders) and the demultiplex version of each pool from `final_demux_versions.tsv` instead of `demultiplex_paths`. Today it reads the combined `_final_info.tsv`.
- [ ] `final_demux_versions.tsv` for multi-VCF runs (needs the `vcf_type` of each pool; see Compatibility → Multi-VCF).

---

## Demultiplexing

- [x] Option to reuse a previous `cellSNP` run for `vireoSNP` (`gt_demux_pipeline.reuse_cellsnp_dir`); useful when there's a vcf for a w_gt run but the SNPs haven't changed (many WGS).
- [x] Option to use a per-pool vcf only for pileups and not for vireo (`gt_demux_pipeline.pool_vcf_use`).
- [x] Retain variant matching stats (variants matched to donor VCF out of total pileup): `vireoSNP` writes `vireo.log`, `demul_samples.py` adds the counts to the demux info and the MAP files map them for `update_logs.py`.
- [x] Beautify `get_filt_barcodes` (now in `lib/io.py`).

---

## Wildcards & modality (one refactor)

- [x] `wildcards.modality` handles cDNA/RNA or ATAC rather than functions: the bams of each modality are looked up in `lib/modality.py`.
- [x] Simplify the wildcard structure: global `wildcard_constraints` in the `Snakefile` (`pool`, `donor`, `modality`, and `vcf_type` limited to the vcf columns of `vcf_info`). The `fold_struct_*` keys stay as they are until the config is restructured (Needs your input).
- [x] Handle varying output dirs for:
  - Demux output (single method vs. simultaneous both): see Demultiplexing outputs & stats
  - BAM splitting (finalizing vs. genotype purposes): `gt_check/` and `final/` sub-folders
- [x] Use the `branch` function to handle rules with varied inputs (`create_inp_cellSNP`, `create_inp_splitBams`, `split_bams`).
- [x] Simplify target functions in `produce_targets`: a table of stages per `last_step`; unsupported combinations (multiome with hashsolo, multi-VCF beyond vireo or with split bams) raise a clear error.

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
  - [ ] `split_bams.smk` (entire split bam pipeline): only the cDNA bam is split for now (`bamfilt_by_CB` and `split_bams` have no `{modality}`)
  - [x] `filt_chr_bams` (→ `filt_chr_bams_multiome_cdna` / `_atac`, one loop over the modalities)
- [x] Permanent fix for rule `filt_chr_bams_multiome`.
- [ ] Multi-VCF: cellSNP's regions and vireo's donor genotypes don't vary by `vcf_type` yet (only the paths do); hashsolo and split bams aren't supported.
- [ ] Multiome with hashsolo (`solo` / `both` / `add_*`).
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

## Needs your input

These can't be done from the code alone; each needs a decision or details from you.

- [ ] Restructure `new_config.yaml`: pick the new layout. It renames most keys, so every existing config must be migrated (a migration script could convert old configs). Separate:
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
- [ ] Move pool calculation logic from `input_processing.smk` into config, inspired by `downloadFromSynapse.py`: that script isn't in this repo, so share it or describe the logic.
- [ ] Handle pools whose `cellSNP.base.vcf.gz` is empty: which failure do you see? vireo fails before `create_inp_splitBams` runs, so either make cellSNP/vireo a checkpoint or have vireo write a sentinel the checkpoint reads.
- [ ] Fix `demultiplex_helper_funcs.py` for double HTOs in `parse_file`: what input fails, and how?
- [ ] Resolve the multi-HTO setup in `kite.smk`: what currently fails?

---

## Benched

- [ ] Dash app: pair-wise sample comparison (donor vs donor or pool vs pool) sub-plot pop up.
- [ ] Agentic AI support for genotype-based demultiplexing using Ollama (with support for public vLLMs).
- [ ] Shannon entropy for a better selection of SNPs (benchmark against the existing `max50k` / `25pct` sets first).
