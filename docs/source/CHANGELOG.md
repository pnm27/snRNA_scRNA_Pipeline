# Changelog

## Pipeline & Workflow

- Renamed workflows and sub-workflows for consistency:
  - `pheno_demux3.snkmk` → `genotype_demux.snkmk`
  - `demultiplex_no_argp.snkmk` → `demultiplex.snkmk`
  - `demul_samples_no_argp.py` → `demul_samples.py`
- Removed unused sub-snakefiles: `helper_functions.smk` and `analyse_vireo.smk`
- Removed dependence on `resources.smk`; all resource requirements are now defined within each Snakemake file
- Added support for Snakemake ≥ v8:
  - Added `workflow_profile/config.yaml`, reflected in `run_snakemake.sh`
  - Replaced `threads` directive with `resources: cpus_per_task`
  - `lsf.yaml` retained for Snakemake < v8
- Standardized single wildcard name to `pool` (previously mixed: `num`, `id1`, `id2`); double wildcards retained
- Added `multiome_alignment` as a new module with `cellranger.smk` (currently supports cellranger arc count only)

---

## Demultiplexing (`demultiplex.smk` / `demul_samples.py`)

### Rule restructuring
- Split the single `demux_samples` rule into three rules reflecting output type:
  - `demux_samples_solo`, `demux_samples_vireo`, `demux_samples_both`
- Added `demultiplex` config section for `demux_samples_both`
- Added `create_h5ad_only` rule for experiments not requiring demultiplexing (no new config option needed)
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
- Changed `vcf_type` wildcard to support both multi-vcf and multiome setups
- Function `ret_htos_calico_solo` now returns a tuple instead of a list
- hashsolo now runs via `scanpy.external` (deprecated: `solo-sc`)

---

## Genotype Demultiplexing (`genotype_demux.smk`)

### cellSNP
- Simplified rule; now only runs `cellsnp-lite` (SNP filtering is better handled outside the pipeline)
- Added option to run without a reference VCF (1000 Genomes Project VCF is the minimum requirement)
- Added option to skip subsetting reference SNP VCFs
- Fixed UMI tag selection for multiome-ATAC
- Fixed barcode handling for multiome (appends `-1` suffix to match cellranger BAM barcodes)
- Replaced indexed array shell logic with `cmd_str_csnp` parameter function (was broken in Snakemake ≥ v8)

### vireoSNP
- Replaced indexed array shell logic with `cmd_str_vireo` parameter function (was broken in Snakemake ≥ v8)
- `vcf_info` file is now expected to have headers
- Simplified command-line execution

### `vcf_info` file specification
- **Must have headers**
- Accepted column orders:
  - Without genotypes: `pool, n_donors`
  - With genotypes (partial or full): `pool, n_donors, donor_names, vcf[, vcf2, vcf3, ...]`
- Deprecated `vcf_info_columns`; column names are free but order is fixed as above
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
- `targets_multiome` now uses a dispatch dict instead of an if/elif chain
- Both dispatcher functions now raise on unrecognised steps (previously returned empty lists silently)
- Moved the following to a new function to accommodate checkpointed `split_bams`:
  - `starsolo_split_bams`, `starsolo_split_bams_gt_demux`, `starsolo_split_bams_gt_demux_multi_vcf`, `starsolo_gt_demux_identify_swaps`
- Fixed logic issue for single-cell/nucleus processing
- Removed dead code and last_step branches: `targets_resolve_swaps_gt_demux`, `targets_multibamsummary`, `targets_multibamsummaryPlotCorr`, `targets_all`

---

## Logging & QC (`update_logs.py` / `run_update_logs.sh`)

### `update_logs.py`
- All parameters except `map_file`, `output_file`, and `bam_dir` are now optional (default to `None`)
- Missing `picard_dir` → skips GCBias and RNASeq Metrics collection
- Missing `demul_dir` → skips demultiplexing info collection
- Added annotation support via the same JSON file used in `demul_samples.py`
- Refactored internals: 1 loop for conversion, 1 loop for calculations, fully vectorized with `pd.to_numeric`
- No duplicated code; structured outputs; easier to extend
- `pd.concat` usage updated per FutureWarning
- Backward support partially removed (advised to create new output files)
- Added support for STARsolo 2.7.10 via `Final_out_MAP_2_7_10a_latest.tsv`
- New reference file: `Final_out_MAP_2_7_10a_latest_info.xlsx`

### `run_update_logs.sh`
- Parameters and arguments now arranged as associative arrays
- Accepts new parameters for annotations: wet lab file and annotation JSON
- Use `empty` values to emulate missing picard or demultiplexing data

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
- Deprecated `vcf_info_columns`

---

## STARsolo

- Per-pool memory allocation for `STARsolo_sort` based on input type (`wildcards.pool`)
- Modelled resource consumption approach for STARsolo
- `samtools index` is now a separate rule to prevent deletion of `STARsolo_sort` outputs on indexing errors
- Input file lines are now sorted when accessed via path in `input_processing.smk`

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
