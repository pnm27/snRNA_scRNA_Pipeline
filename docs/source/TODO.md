# TODO

Open items only: completed ones are listed in `CHANGELOG.md`.

---

## Config restructure

New layout (`helper_py_scripts/migrate_config.py` has the full key mapping): `references` (genomes, barcodes, regions, tool versions), `workflow` (aligner, until, picard_metrics, gt_check), `experiment` (inputs, hto, demultiplexing, qc), `modules` (the `*_pipeline` blocks, each with its `fold_struct*`) and `execution` (rule `resources`; site defaults, Minerva first). `last_step` becomes `workflow.aligner` + `workflow.until` (+ `experiment.demultiplexing.genotype.multi_vcf`); `demux_run_type` is absorbed by `run_tag`, with `experiment.demultiplexing.genotype.enabled` replacing the `wo_gt` checks in `lib/io.py`. Old configs are converted by the migration script only. Stages 1–3 are done (see CHANGELOG); the caveats are in `config_caveats.md`.

- [ ] The shared `fold_struct` (`Sample_{pool}-cDNA/{pool}-cDNA`) moved to `modules.STARsolo.fold_struct`, but other modules still use it: `kite.smk` (the fastq folder, with `-cDNA` → `-HTO`), `lib/io.py` and `lib/params.py` (fastq folders), `lib/config_models.py` (`modules.gt_demux`'s h5ad without demultiplexing). Give each of them its own `fold_struct` so only STARsolo uses `modules.STARsolo.fold_struct`.

---

## Testing

- [ ] Add a tiny dataset to `tests/scRNAseq` and `tests/pooled_snRNAseq` (both empty) for real, not only dry, runs.
- [ ] Run the dry-run and unit tests in GitHub Actions.

---

## Demultiplexing outputs & stats

The layout is described in `docs/source/sub_workflows.md` (Demultiplexing outputs).

- [ ] `update_logs.py`: read the basic/hashsolo/vireo stats files separately (its `Basic_Demultiplex_stats` etc. placeholders) and the demultiplex version of each pool from `final_demux_versions.tsv` instead of `demultiplex_paths`. Today it reads the combined `_final_info.tsv`.

---

## Aligners

- [ ] hashsolo with `cellranger count`: its barcodes end with `-1` and the HTO counts' (kite) don't, so they need matching (raises an error for now).
- [ ] Picard metrics for `cellranger count` and multiome runs: Picard runs on STARsolo runs only (it uses STARsolo's `genome_pick` reference).
- [ ] `update_logs.py`: compile `cellranger count`'s `metrics_summary.csv`.

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
- [ ] `update_logs.py` needs `--common_annotations`: without it, `data`, `cols_list` and the demux paths are never set and it crashes. Make the annotations optional (and read the annotation schema relative to the script, not the working directory).
- [ ] Split `demultiplex_helper_funcs.py` (~1,600 lines) into smaller modules, after the config restructure.

---

## Experimental / Future

- [ ] Add new Picard metrics
- [ ] Incorporate git submodules for external repos (e.g. kite `featuremap`)

---

## Documentation

- [ ] `sub_workflows.md`: the Module description table still lists the old `last_step` values (and its rows are placeholders); describe the stages of `workflow: until` instead.
- [ ] Add tutorials (they wait on items in Needs your input):
  - Pooled snRNA-seq and scRNA-seq with multiple wildcards (the YAML input of `input_processing.smk` is still being reworked)
  - Double HTOs (see the double- and multi-HTO items)

---

## Needs your input

These can't be done from the code alone; each needs a decision or details from you.

- [ ] Move pool calculation logic from `input_processing.smk` into config, inspired by `downloadFromSynapse.py`: that script isn't in this repo, so share it or describe the logic.
- [ ] Handle pools whose `cellSNP.base.vcf.gz` is empty: which failure do you see? vireo fails before `create_inp_splitBams` runs, so either make cellSNP/vireo a checkpoint or have vireo write a sentinel the checkpoint reads.
- [ ] Fix `demultiplex_helper_funcs.py` for double HTOs in `parse_file`: what input fails, and how?
- [ ] Resolve the multi-HTO setup in `kite.smk`: what currently fails?

---

## Benched

- [ ] Dash app: pair-wise sample comparison (donor vs donor or pool vs pool) sub-plot pop up.
- [ ] Agentic AI support for genotype-based demultiplexing using Ollama (with support for public vLLMs).
- [ ] Shannon entropy for a better selection of SNPs (benchmark against the existing `max50k` / `25pct` sets first).
