# TODO

## Architecture & Design

- [ ] Currently, the creation of h5ad from vireo output lies within the rule `demux_samples`. Split that rule into 2 rules:
  - [ ] RULE A: An h5ad is created after vireo output
  - [ ] RULE B: `demux_samples` picks proper version of vireo demultiplexing, corrects swaps and creates a final h5ad.
- [ ] Currently, the dash app can't compare samples. Make it so that a pair-wise sample (can be donor vs donor or pool vs pool) sub-plot pop up.
- [ ] Next update will be split the design for sc/snRNA seq, scATAC and multiome. Pipeline will be restructured.
- [ ] `wildcards.modality` will handle cDNA/RNA, HTO or ATAC rather than functions (Input, params, etc.)
- [ ] Split the demultiplexing info. It should be like basic stats (common file), hashsolo to hashsolo_stats file and vireo to vireo_stats file. This would make all things distinct (stats saved per demux per config/params). Not only will it be organized well but also mark the final version for each pool (by specifying in swap_correction_df).
- [ ] Use `branch` function to handle rules with varied inputs
- [ ] Move pool calculation logic from `input_processing.smk` into config (inspired by `downloadFromSynapse.py`)
- [ ] Change how demultiplexing directory variants are handled — replace separate config dirs (e.g. `wo_gt`, `w_gt`) with a suffix set at the top of the config file
- [ ] Employ a clear strategy for `final_count_matrix_dir` across demux scenarios:
  - Both tools run simultaneously
  - Tools run sequentially (name each run separately or record the order)
- [ ] Simplify wildcard structure so folder structure encodes wildcards naturally
- [ ] Revamp wildcards to correctly handle varying output dirs for:
  - Demux output (single method vs. simultaneous both)
  - BAM splitting (finalizing vs. genotype purposes)
- [ ] Remove dependency on STARsolo as the only aligner
- [ ] Add input function for rules using `genefull_matrices` to select either `Gene` or `GeneFull` per project
- [ ] Absorb HTO, cDNA and ATAC into the wildcard `modality`.
- [ ] Restructure `new_config.yaml` to:
  - separate:
    - [ ]  static resources
    - [ ]  experiment definitions
    - [ ]  workflow behavior:
      **OLD STYLE**

      ```sh
      demux_type: solo
      demux_run_type: wo_gt
      hto_demux_type: single
      ```

      **NEW STYLE**

      ```sh
      experiment:
        demultiplexing:
          method: solo
          genotype:
            enabled: false

      hto:
        enabled: true
        multiplexing: single
      ```

    - [ ]  execution modules
    - [ ]  Output naming conventions
- [ ]  Support ArchR and Signac QC for ATAC QC metrics.

---

## Compatibility

- [ ] Make the pipeline compatible with:
  - [ ] Multiome
  - [ ] Multi-module
  - [ ] Multi-VCF
- [ ] Generalize outputs to accommodate multiome:
  - [ ] `split_bams.smk` (entire split bam pipeline)
  - [ ] `filt_chr_bams` (→ `filt_chr_bams_multiome`)
- [ ] Add `CollectInsertSizeMetrics` for the ATAC modality in multiome (`picard_metrics`)
- [ ] Add multiome support in `update_logs.py`
- [ ] Permanent fix for rule `filt_chr_bams_multiome`

---

## Demultiplexing

- [ ] Fix `vcf_type` wildcard issue in the rule `demux_samples`
- [ ] Fix automatic output directory selection in `demux_samples`
- [ ] When adding calico_solo or vireo results, include the existing demux stats file as input and append to it
- [ ] Retain vireo rerun information in both the demultiplex info file and `update_logs`:
  - [x] *swap correction* addresses proper demultiplex file per each pool.
  - [ ] *swap correction* has appropriately identified donors per each pool.
- [ ] Retain variant matching stats from log runs (variants matched to donor VCF out of total pileup):

  ```shell
  find -mindepth 1 -maxdepth 1 -type d -exec sh -c \
    'a=$(sed "s#\./##g" <<< {}); b=$(ls -ltr ${a} | tail -1 | rev | cut -d " " -f1 | rev); \
     grep "variants matched to donor VCF" ${a}/${b}' \;
  ```

- [ ] Fix `demultiplex_helper_funcs.py` for double HTOs in `parse_file`
- [ ] Beautify `get_filt_barcodes` in `genotype_demux.smk`

---

## cellranger

- [ ] Add support for cellranger.

---

## STARsolo

- [ ] Add WASP mode

---

## Code Quality

- [ ] Add Python logging support using the `logging` module
- [ ] Simplify target functions in `produce_targets`
- [x] Remove `run` directive from `kite`, `pheno_demux3`, and `STARsolo`
- [x] Combine sub-workflows `split_bams` and `split_bams_gt`
- [x] Add config options for extra software parameters

---

## Experimental / Future

- [ ] Agentic AI support for genotype-based demultiplexing using Ollama (with support for public vLLMs)
- [ ] Add new Picard metrics
- [ ] Incorporate git submodules for external repos

---

## Documentation & Schemas

- [ ] Write schemas
- [ ] Add tutorials:
  - Pooled snRNA-seq: single wildcard, multiple wildcards
  - scRNA-seq: single wildcard, multiple wildcards
  - Double HTOs
- [ ] Search ranking config for ReadTheDocs
