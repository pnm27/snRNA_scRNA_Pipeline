# TODO

Open items only: completed ones are listed in `CHANGELOG.md`.

---

## Groundwork: config & schemas

- [ ] Tighten the schema (`schema/config.schema.yaml`) once the config is restructured (required nested keys, enums for `last_step`).
- [ ] Consider growing `lib/config_models.py` (`DemuxConfig`) into typed config.

---

## Demultiplexing outputs & stats

The layout is described in `docs/source/sub_workflows.md` (Demultiplexing outputs).

- [ ] `update_logs.py`: read the basic/hashsolo/vireo stats files separately (its `Basic_Demultiplex_stats` etc. placeholders) and the demultiplex version of each pool from `final_demux_versions.tsv` instead of `demultiplex_paths`. Today it reads the combined `_final_info.tsv`.

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
- [ ] Generalize `split_bams.smk` (the entire split bam pipeline) for multiome: only the cDNA bam is split for now (`bamfilt_by_CB` and `split_bams` have no `{modality}`).
- [ ] Multi-VCF: cellSNP's regions and vireo's donor genotypes don't vary by `vcf_type` yet (only the paths do); hashsolo and split bams aren't supported.
- [ ] Multiome with hashsolo (`solo` / `both` / `add_*`).
- [ ] Add `CollectInsertSizeMetrics` for the ATAC modality in multiome (`picard_metrics`).
- [ ] Add multiome support in `update_logs.py`.
- [ ] Support ArchR and Signac QC for ATAC QC metrics.

---

## Code Quality

- [ ] Python logging with the `logging` module in the remaining helper scripts (`update_logs.py` already uses it).

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
