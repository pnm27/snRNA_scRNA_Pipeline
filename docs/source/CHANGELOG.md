# Changelog

## Pipeline & Workflow

- Renamed workflows and sub-workflows for consistency:
  - `pheno_demux3.snkmk` → `genotype_demux.snkmk`
  - `demultiplex_no_argp.snkmk` → `demultiplex.snkmk`
  - `demul_samples_no_argp.py` → `demul_samples.py`
- Removed unused sub-snakefiles: `helper_functions.smk` and `analyse_vireo.smk`
- All resource requirements are now defined in one table in `rules/resources.smk`; rules use `**get_resources("<rule>")` (int, or `base + step * (attempt - 1)` with an optional `max`), and the config's `resources` block can override any rule. `cellSNP` / `vireoSNP` values can differ by VCF suffix (`by_vcf_suffix`). Removed the duplicated and unused `allocate_*` functions and the `# threads` comments for Snakemake < v8
- The config is validated against `schema/config.schema.yaml` at the start of every run; defaults fill in optional keys missing from older configs (e.g. `swap_correction_df`)
- Added dry-run tests (`tests/dryrun/run_dryruns.sh`) that check the DAG of each scenario in `tests/dryrun/scenarios/` without data or a cluster
- Global `wildcard_constraints` in the `Snakefile`: `pool` and `donor` can't contain `/`, `modality` is `cDNA|ATAC|HTO`, and `vcf_type` is one of the vcf columns of `vcf_info`
- Per-modality bam paths are looked up in `lib/modality.py` (`get_bam`, `get_cellsnp_inputs`, `get_umiTag`); the `filt_chr_bams*` rules are one loop over the modalities with a bam
- `create_inp_cellSNP`, `create_inp_splitBams` and `split_bams` pick their input with `branch()` (replacing `get_filt_barcodes`, `get_inp_splitBam` and `get_bam_to_split`)
- `rule all`'s targets come from a table of stages per `last_step` (`rules/produce_targets.smk`); combinations not supported yet raise a clear error
- Split bams (donor hash file, per-donor bams, logs, temp dir) go to `gt_check/` or `final/` sub-folders; `split_bams_dir2`, `per_donor_split_log_dir2`, `new_temp_dir2` and the `split_bams_proxy_dir*` keys are no longer used. Per-donor bams from earlier runs are outside these sub-folders, so move them in to avoid re-creating them
- Multiome split bams use vireo's cDNA output and the cDNA bam
- Fixed: `STARsolo_cellsnp` and split bams targets pointing at files no rule makes; multi-VCF targets (undefined `VCF_TYPE`); `identify_swaps` running only for the first donor of each pool; missing `SimpleNamespace` import for multiome; `bamfilt_by_CB` sharing one `{pool}_bc.txt` across pools
- Dry-run tests can mark the jobs up to the donor checkpoint as done (`# donors:`) to check per-donor jobs, and check expected errors (`# fails:`)
- Unit test for STARsolo's run time (`tests/unit/`)
- Added support for Snakemake ≥ v8:
  - Added `workflow_profile/config.yaml`, reflected in `run_snakemake.sh`
  - Replaced `threads` directive with `resources: cpus_per_task`
  - `lsf.yaml` retained for Snakemake < v8
- Standardized single wildcard name to `pool` (previously mixed: `num`, `id1`, `id2`); double wildcards retained
- Added `multiome_alignment` as a new module with `cellranger.smk` (currently supports cellranger arc count only)
- Consolidate demultiplexing into a single rule (combine all 4 `demux_samples_*` rules)
- Removed the `run` directives from the `kite`, `genotype_demux` (was `pheno_demux3`) and `STARsolo` rules
- Config options for extra command-line parameters of the tools (`extra_params`, `mbv_extra_opt`)
- Most of the parameter functions have been moved to `lib/params.py`.
- Most of the input functions have been moved to `lib/io.py`.
- Input file now supports a 2-columned (tab-separated) file containing *pool*, *SS_MEM*, where *SS_MEM* is RAM memory per thread/task.

---

## Demultiplexing (`demultiplex.smk` / `demul_samples.py`)

### Rule restructuring

- Split `demux_samples` into `demux_run_<tool>` (one h5ad per tool: `solo`, `vireo`, or `none` without demultiplexing) and `demux_final` (the final h5ad of each pool), with output paths in `lib/config_models.py` (`_output_cfg` removed)
  - `demux_final` combines the runs in the pool's row of the new swap correction file (`demultiplex.swap_correction`, versions as `<tool>:<run_tag>`), else this run's tools; it adds a `Donor_final` column, renamed by the swap correction file
  - `add_solo` / `add_vireo` take the other tool's run from `demultiplex.previous_run_tag` and no longer start from the earlier final h5ad
  - `demultiplex.demux_count_matrix_dir` (was only for `both`) now holds the final h5ads of all runs; new `fold_struct_final` (default `{pool}`)
  - `gt_demux_pipeline.swap_correction_df` is no longer used
  - With `create_final_h5ad: False`, solo runs now also stop at calico_solo's output, as vireo runs did
  - Resources: `demux_samples` → `demux_run` and `demux_final` (update `resources:` overrides in configs)
- New `run_tag` config key: each demultiplexing run writes to its own sub-folder of the calico_solo, cellSNP, vireoSNP, final count matrix and demultiplex info dirs
- Split the demultiplex info: `demux_run_<tool>` writes only its tool's stats (`_hashsolo_stats.tsv` / `_vireo_stats.tsv`) and keeps the QC stats in the h5ad (`uns['qc_stats']`); `demux_final` writes `_basic_stats.tsv` (QC and the versions used) and the combined info read by `update_logs.py`
- Swap correction versions can name runs of two tools, e.g. `solo:hto,vireo:w_gt`
- New rule `collect_demux_versions` writes `final_demux_versions.tsv` (demultiplex versions, swap correction, donors and cells passing QC per pool); multi-VCF runs get one table per `vcf_type` (`final_demux_versions_<vcf_type>.tsv`)
- Multiome demultiplexing runs once per pool (`demux_samples` ran once per modality and wrote the same h5ad twice)
- `split_bams_pipeline.split_by.input: h5ad` reads the h5ad of `split_by.demux` (it always read hashsolo's); new `input: final` splits by the final h5ad (`Donor_final`)
- Fixed `demul_samples.py` with calico_solo: `cols` used before being set, the classification run twice, and `args.pref` (doesn't exist)
- Fixed `create_inp_splitBam_consolidated.py` with an h5ad input crashing without `--demux_suffix` and keeping hashsolo's doublets and negatives
- Split the single `demux_samples` rule into three rules reflecting output type:
  - `demux_samples_solo`, `demux_samples_vireo`, `demux_samples_both`
- Added `demultiplex` config section for `demux_samples_both`
- Added `create_h5ad_only` rule for experiments not requiring demultiplexing (no new config option needed)
- Fixed `demux_type: null` (h5ad without demultiplexing) crashing; null, `none` and empty values are now handled by `get_demux_type()` in `lib/utils.py`
- Fixed `demux_type: both` targets not matching the single combined h5ad written by `demux_samples`; its demux info file now goes to `demultiplex.demultiplex_info_dir`
- Fixed doubled `{vcf_type}` directories in multi-VCF runs; `demux_samples` outputs now carry `_{vcf_type}` when their fold struct lacks it
- Fixed multiome `demux_samples` passing the same modality's vireo output twice (now ATAC and cDNA)
- Made `demux_info` parameter optional (previously positional) in the rule handling addition of new demux to a final count matrix
- Simplified shell scripts in `demux_samples` and `add_obs_to_final_count_matrix`
- `demux_info` parameter changed from positional to optional in the add-demux rule

### Functional additions

- Added `add_solo` / `add_vireo` config options to support running a second demultiplexing tool on top of an existing result
- Added config option to enable/disable h5ad creation during demultiplexing (useful as a switch during genotype checks)
- `demul_samples.py` now saves QC-failed genes in the final h5ad under `var` column `QC_pass`
- Added multiome demultiplexing support in `demux_samples_vireo` (set global `ONLY_VIREO`)
- Added multiome support in `demul_samples.py` via append mode for vireo
- `demultiplex_helper_funcs.py`:
  - Added support for annotating outputs using a JSON file
  - Created dataclass for custom outputs and views
  - Added a new parameter to the function `get_filename` to support latest demultiplexing info filenames.
  - Changed previously named function `process_swap_correction` to `final_demultiplex_version` due to the name being a misnomer.
  - Added a new `process_swap_correction` function which now processes donor name correction.
- Changed `vcf_type` wildcard to support both multi-vcf and multiome setups
- Function `ret_htos_calico_solo` now returns a tuple instead of a list
- `demul_samples.py` adds vireo's variant counts (pileup, matched to donor VCF, used) to the demux info, read from `vireo.log` next to `donor_ids.tsv` (function `vireo_log_stats`)
- hashsolo now runs via `scanpy.external` (deprecated: `solo-sc`)
- Add `--keep-going` parameter in `run_snakemake.sh` so that independent jobs get submitted even if some other job fails.

---

## Genotype Demultiplexing (`genotype_demux.smk`)

### cellSNP

- Simplified rule; now only runs `cellsnp-lite` (SNP filtering is better handled outside the pipeline)
- Added option to run without a reference VCF (1000 Genomes Project VCF is the minimum requirement)
- Added option to skip subsetting reference SNP VCFs
- Fixed UMI tag selection for multiome-ATAC
- Fixed barcode handling for multiome (appends `-1` suffix to match cellranger BAM barcodes)
- Replaced indexed array shell logic with `cmd_str_csnp` parameter function (was broken in Snakemake ≥ v8)
- Added `gt_demux_pipeline.pool_vcf_use` to use a per-pool vcf only for the pileup and/or for vireo (default: both)
- Memory/time scaling by `vcf_type` suffix now applies (`'vcf_type' in wildcards` checked wildcard values, so it never did)

### vireoSNP

- Replaced indexed array shell logic with `cmd_str_vireo` parameter function (was broken in Snakemake ≥ v8)
- `vcf_info` file is now expected to have headers
- Simplified command-line execution
- Added `gt_demux_pipeline.reuse_cellsnp_dir` to run vireo on the cellSNP outputs of a previous run
- vireo's output is now also written to `vireo.log` (rule `log`), next to `donor_ids.tsv`
- With a 2-column `vcf_info` and a genotyped run, vireo now gets `donors_per_pool_file` as donor genotypes (was the 1kGP sites VCF)
- Clear error when `vcf_info` (needed for `-N`) is missing

### `vcf_info` file specification

- **Must have headers**
- Accepted column orders:
  - Without genotypes: `pool, n_donors`
  - With genotypes (partial or full): `pool, n_donors, donor_names, vcf[, vcf2, vcf3, ...]`
- Deprecated `vcf_info_columns`, except `vcf_info_columns: vcf`, which names the vcf columns of multi-VCF runs (the values of `vcf_type`); column names are free but order is fixed as above
- `precursor_to_perPoolVcf.ipynb` now templates the `vcf_info` file

---

## BAM Splitting (`split_bams.smk`)

- Consolidated GT and non-GT versions of `split_bams` script
- `mito_file` moved from `output` to `params`; `bed_file` moved from `input` to `params`
- Major revisions to `create_per_donor_bams.bash`:
  - Consolidated GT and non-GT versions
  - More elegant mito file handling
  - Added argument parsing (with backward compatibility for positional args)
  - Input directories no longer need to follow specific logic, but must end with `/`
- `create_inp_splitBams` rule:
  - Uses consolidated script for barcode file creation from both h5ad and raw files
  - Removed overwrite option from the rule (still present in the Python script)
  - Supports single demux (script can handle multiple)
  - h5ad-only input supported for calico_solo; vireo output supported as-is
- Added multiome BAM splitting support (cDNA and ATAC modalities)
- Fixed wrong inputs in `split_bams.smk`
- Fixed `get_inp_splitBam` returning nothing when splitting by raw hashsolo output; unknown `split_by.demux` values now raise an error
- Added wildcard constraint for `qtltools_mbv` to handle:
  - Multiome: `<path>/BD2-Set-4-a/cDNA_Pitt-DNA-PFC-694.bamstat.txt`
  - Standard: `BD2-Set-10-a/donor4.bamstat.txt`

---

## Target Production (`produce_targets.smk`)

- Eliminated `wc_d`; all `expand(..., zip, **wc_d)` calls replaced with:
  - `_expand_pools(template)` — static, pool-only expansion
  - `_expand_pools_donors(template, wildcards)` — runtime, pool+donor expansion
- Target builder functions now consistently return lists
- Extracted `_picard_targets` to eliminate duplicated config key lookups
- Both dispatcher functions now raise on unrecognised steps (previously returned empty lists silently)
- Moved the following to a new function to accommodate checkpointed `split_bams`:
  - `starsolo_split_bams`, `starsolo_split_bams_gt_demux`, `starsolo_split_bams_gt_demux_multi_vcf`, `starsolo_gt_demux_identify_swaps`
- Fixed logic issue for single-cell/nucleus processing
- Removed dead code and last_step branches: `targets_resolve_swaps_gt_demux`, `targets_multibamsummary`, `targets_multibamsummaryPlotCorr`, `targets_all`
- `targets_demux` takes no arguments; its targets come from `DemuxConfig` (`lib/config_models.py`)
- Fixed missing `modality` expansion in the `STARsolo_kb_solo` target
  
---

## Logging & QC (`update_logs.py` / `run_update_logs.sh`)

### `update_logs.py`

- All parameters except `map_file`, `output_file`, and `bam_dir` are now optional (default to `None`)
- Missing `picard_dir` → skips GCBias and RNASeq Metrics collection
- Missing `demul_dir` → skips demultiplexing info collection
- Added annotation support via the same JSON file used in `demul_samples.py`
- Refactored internals: 1 loop for conversion, 1 loop for calculations, fully vectorized with `pd.to_numeric`.
- No duplicated code; structured outputs; easier to extend.
- Function `write_logs` moved to `demultiplex_helper_funcs.py`.
- Logging is now supported and can be modified using *verbose* argument.
- *Swap correction* is now supported. Metrics from the appropriate run of each pool is used for compilation.
- `pd.concat` usage updated per FutureWarning.
- Backward support partially removed (advised to create new output files).
- Added support for STARsolo 2.7.10 via `Final_out_MAP_2_7_10a_latest.tsv`.
- New reference file: `Final_out_MAP_2_7_10a_latest_info.xlsx`.
- Now, when STARsolo's log file is not present it will still collect info while providing a warning.
- Added verbosity as an option to `run_update_logs.sh`.

### `run_update_logs.sh`

- Parameters and arguments now arranged as associative arrays.
- Accepts new parameters for annotations: wet lab file and annotation JSON.
- Use `empty` values to emulate missing picard or demultiplexing data.
- Supports *swap correction*.
- To support newer demultiplexing info output names (*cDNA_STARsolo_vS_info.tsv*, *ATAC_STARsolo_vS_info.tsv*).
- Changed the function `process_swap_correction`
- Builds the command as a bash array (paths with spaces work) and passes samples with `--samples`
- `check_required_args` fails fast when a required arg or `INPUT_FILE` is missing
- Output suffixes are set in one `opt_values` map (empty = default suffix, remove a line to skip that output); unknown keys are warned about and skipped
- MAP files (`Final_out_MAP_*_latest*.tsv`) map the new vireo variant counts: `N_VARIANTS_PILEUP_VS`, `N_VARIANTS_MATCHED_VS`, `N_VARIANTS_USED_VS`
- STARsolo's run time: new `(STAR, REG, RUN_TIME)` column after `JOB_FINISH` (`JOB_FINISH - JOB_START` as HH:MM:SS) in the latest MAP files, computed once all samples are compiled; reads `2025-05-13 19:51:55` and STAR's year-less times

---

## Wet Lab Info (`create_wet_lab_info.py`)

- Can now run without a converter file
- Saves donor file alongside the wet lab compilation file
- argparse fully documented
- Mirrors actions for donor and multiplex compilations
- Fixed file extension issue when reading existing wet lab info files for updates
- Fixed issue where some lines from certain files (e.g. AMP) were missed

---

## Configuration (`new_config.yaml`)

- Renamed config keys:
  - `demux_pipeline` → `hashsolo_demux_pipeline`
  - `split_bams_pipeline_gt_demux` → `split_bams_pipeline_gt`
  - `gt_conv` → `file` (under `donorName_conv` in `gt_demux_pipeline`)
  - `mito` → `mito_prefix`
  - `gtf_file` → `gtf`, `genome_fasta` → `fasta`, `genome_dir` → `genome`, `sjdboverhang` → `overhang` (all resolved via `genome_pick`)
- Added `demultiplex` section for `demux_samples_both`
- Added PICARD option
- Fixed genome/GTF/FASTA incompatibility: selecting `genome_pick` in `STARsolo_pipeline` now automatically resolves `gtf`, `fasta`, `genome`, `overhang`, and `gene_info_file` via anchors and references
- Deprecated `vcf_info_columns` (except `vcf_info_columns: vcf`, for multi-VCF runs).
- Removed `cellranger_arc_ref` and replaced it with genome selection (similar to STARsolo's)
- Added `cellranger_arc_count.ref` to `new_config.example.yaml` (it failed to load without it)
- New optional keys: `gt_demux_pipeline.reuse_cellsnp_dir`, `gt_demux_pipeline.pool_vcf_use` and top-level `resources`
- Config restructure, stage 1: `helper_py_scripts/migrate_config.py` converts configs to the new layout (`references`, `workflow`, `experiment`, `modules`, `execution`; see TODO). It keeps comments and anchors and also converts partial configs such as the dry-run scenarios. `last_step` becomes `workflow.aligner` + `workflow.until` (+ `genotype.multi_vcf`), `demux_type` becomes `demultiplexing.method` + `add_to_previous`, `demux_run_type` becomes `genotype.enabled`, `hto_demux_type` becomes `hto.enabled` + `multiplexing`, the `fold_struct_*` keys move into their modules (`fold_struct` into `modules.STARsolo`), and the unused `fold_struct_gt_demux2` is dropped. It stops on duplicate keys, since snakemake keeps the last value. The pipeline still reads the old layout until stage 2.
- Unit tests for the migration (`tests/unit/test_migrate_config.py`): no key of `new_config.example.yaml` or of any dry-run scenario is lost, the new keys mean what `last_step`, `demux_type`, `demux_run_type` and `hto_demux_type` did, and migrating a scenario apart from the example gives the same merged config
- Removed the duplicate `STARsolo_pipeline.gene_barcodes` from `new_config.example.yaml` (snakemake used the second, `.tsv.gz`)
- Config restructure, stage 2: the pipeline reads the new layout (`config['modules']['STARsolo'][...]`, etc.), and `new_config.example.yaml`, the dry-run scenarios and the tutorials use it. A config in the old layout stops the run with the command to convert it (dry-run scenario `old_layout`)
  - `workflow: until` replaces `last_step`'s stages (`rules/produce_targets.smk`: `UNTIL_STAGES`). The alignment is always a target, so runs that used to stop short of it (e.g. `STARsolo_demux` without Picard) now also index the bam, and Picard's metrics run whenever `picard_metrics` is set
  - `workflow: aligner` (`lib/aligner.py`) and `genotype: multi_vcf` (`lib/config_models.py`) replace the prefix and suffix of `last_step`; `get_demux_type` (`lib/utils.py`) reads `method` and `add_to_previous`; `is_genotyped` replaces the `demux_run_type == 'wo_gt'` checks (`lib/io.py`); `SINGLE_HTO` / `MULTI_HTO` replace the `hto_demux_type` checks (`kite.smk`); rule resource overrides are read from `execution: resources`
  - `schema/config.schema.yaml` follows the new layout (the old keys' types, defaults and required keys at their new paths; the duplicate `gene_barcodes`, whose `.tsv.gz` default was overridden, is merged). Tightening it is stage 3
  - `migrate_config.py` keeps the comment lines a config starts with (the dry-run scenarios' checks)
  - Dry runs: `JOBS_DIR` writes each scenario's jobs and output paths, to compare them before and after a refactor. All 28 scenarios have the same jobs as before stage 2, except for the indexed bam of the 4 that used `STARsolo_demux` without Picard (and, for `swap_correction`, the alignment it needs)
  - Unit tests check the migration against old-layout copies of the example config and scenarios (`tests/unit/old_layout/`), and that the new example config passes the schema
- Config restructure, stage 3: the config check explains what's wrong (`docs/source/config_caveats.md` lists the caveats, caught or not)
  - `lib/config_check.py` checks the config before snakemake's `validate()` and lists every problem at once, each with why: the `x-explain` of the schema rule it breaks, or the key's `description`. It also stops on duplicate keys in the config files (YAML keeps the last value) and on the old layout
  - Schema: unknown keys stop the run at the top level and in `workflow`, `experiment` and `execution` (`modules` and `references` stay open); enums for `method`, `add_to_previous` and `hto: multiplexing`; rules across keys (`allOf`): `add_to_previous` must be the other tool and needs `previous_run_tag`, multi-VCF runs (vireo only, `{vcf_type}` in the fold_struct, vcf columns, no split bams), split bams need a demultiplexing method, hashsolo needs HTOs and the wet lab file and isn't supported with cellranger or multiome, HTOs need `multiplexing`, cellranger needs `modules: cellranger`. Warning (`x-warnings`): `picard_metrics` is ignored for cellranger and cellranger-arc
  - `RunSettings` (`lib/utils.py`, `run_settings(config)`) reads the `workflow` and `experiment` keys the rules branch on in one place; `get_aligner`, `get_demux_type`, `is_genotyped`, `is_multi_vcf`, `is_multiome`, kite's `SINGLE_HTO` / `MULTI_HTO`, `until` and `run_tag` use it
  - Tests: `tests/unit/test_config_check.py` (a case per rule: the error and its explanation); dry-run scenarios `unknown_key` and `duplicate_key`. The jobs of the other scenarios are the same as after stage 2

---

## Alignment Support

- The aligner's outputs (count matrix, bams, Picard outputs) come from `lib/aligner.py`, so rules don't build STARsolo paths themselves; the prefix of `last_step` picks the aligner (`starsolo_*`, `cellranger_*` or `multiome_*`)
- `bamfilt_by_CB` writes the aligner's cell-filtered bam (multiome splits read cellranger-arc's, which no rule wrote)

### cellranger

- Added `cellranger count` for GEX-only runs (`last_step: cellranger_*`, config section `cellranger_count`: reference, `include_introns`, `extra_params`); its count matrix and bam are used for demultiplexing (vireo or none) and split bams. hashsolo with it raises an error for now (barcodes end with `-1`)
- Supports ATAC-based (along with RNA-based) vireo demultiplexing via cellranger-arc count.
- Supports cellranger-arc count for multiome preprocessing.

### STARsolo

- `STARsolo_pipeline.count_matrix` picks the count matrix used after the alignment: `GeneFull` (default) or `Gene`; Gene's matrix files are gzipped like GeneFull's and the example config's `gene_*` paths end in `.gz`
- WASP mode (`STARsolo_pipeline.wasp`): `--waspOutputMode SAMtag` with the pool's vcf from `vcf_info`, and `vW` in the SAM attributes. STAR reads only the first sample's genotypes, so rule `wasp_vcf` writes every biallelic SNV of the pool's vcf as heterozygous. Pools without a vcf raise an error
- Per-pool memory allocation for `STARsolo_sort` based on input type (`wildcards.pool`)
- Modelled resource consumption approach for STARsolo
- `samtools index` is now a separate rule to prevent deletion of `STARsolo_sort` outputs on indexing errors
- Input file lines are now sorted when accessed via path in `input_processing.smk`

---

## Documentation

- Added tutorials: scRNA-seq and pooled snRNA-seq (one wildcard)
- ReadTheDocs search ranks the tutorials first and the API reference and development notes last (`search.ranking` in `.readthedocs.yaml`)

---

## Miscellaneous

- Demultiplex info file:
  - Renamed param: `Unique genes` → `gene_ids with an associated gene_name`
  - Added new param for stats on gene IDs lacking an associated gene name
- Accepted input styles for the pipeline (no headers for text files):
  - Pools only
  - Pools with per-pool STARsolo memory (empty = use default)
  - Multi-module setup via YAML (see schema)
- Removed `mode='w+'` when creating outputs in `create_Feat_Barc.py`
- Fixed issues with some packages in `basic_sctools.yaml`
- Upgraded `pyyaml` requirement to `6.0.2` in `sphinxdoc_requirements.txt`
