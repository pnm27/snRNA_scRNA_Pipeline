# Config caveats

The config (`new_config.yaml`) is checked at the start of every run against `schema/config.schema.yaml` by `lib/config_check.py`. It lists every problem at once, each with why it's a problem:

```text
2 problem(s) in the config (see docs/source/config_caveats.md):

experiment: demultiplexing: method: 'solo' isn't allowed here
  Why: hashsolo (method solo or both) isn't supported with cellranger count yet: its cell barcodes end with '-1' and the HTO counts' (kite) don't

Unknown key(s) 'untill' in workflow.
  Why: only these keys are read there: aligner, gt_check, picard_metrics, until. A typo, or a key of the old layout? (helper_py_scripts/migrate_config.py converts old configs)
```

The rules across keys are the `allOf` list of the schema, each with its explanation (`x-explain`); add a rule there to catch a new mistake, and a case to `tests/unit/test_config_check.py`. Warnings (`x-warnings`) are printed and the run goes on.

## Caught at the start of a run

| Caveat | What to do |
|---|---|
| The pipeline reads only the new layout (`references`, `workflow`, `experiment`, `modules`, `execution`); a config with `last_step` stops the run. | `python3 helper_py_scripts/migrate_config.py new_config.yaml -o new_config.migrated.yaml` |
| A key given twice in a YAML file: YAML keeps the last value and silently drops the first (e.g. two `gencode38o90` genomes). The check and `migrate_config.py` both stop on it. Overriding a key merged from an anchor (`<<: *gen38`) is fine. | Keep one of them. |
| Unknown keys stop the run at the top level and in `workflow`, `experiment` (and its sections) and `execution`: a typo would otherwise be ignored. | Fix the key; the error lists the keys read there. |
| `add_to_previous` names the tool of the earlier run, `method` the tool run now: adding hashsolo to a vireo run is `method: solo`, `add_to_previous: vireo` (the old `add_solo`). | Set both, and `modules: demultiplex: previous_run_tag` (the earlier run's `run_tag`; null isn't accepted). |
| Multi-VCF runs (`genotype: multi_vcf: true`) support only `method: vireo`, need `{vcf_type}` in `modules: gt_demux: fold_struct` and the vcf columns in `vcf_info_columns: vcf`, and can't split bams yet. | |
| `until: split_bams` or `identify_swaps` needs a demultiplexing method: bams are split by the donors it finds. | |
| hashsolo (`method: solo` or `both`) needs `experiment: hto: enabled: true` with `multiplexing`, and `experiment: inputs: wet_lab_info`. It isn't supported with cellranger count (barcodes end with `-1`) or multiome (`cellranger_arc`) yet. | |
| With `experiment: hto: enabled: true`, `multiplexing` (`single` or `multiple`) picks the kite rules; without it no HTO rule exists. | |
| `aligner: cellranger` needs `modules: cellranger` (`bams_dir`, `ref`). | |
| Enums: `workflow: aligner` (`STARsolo`, `cellranger`, `cellranger_arc`), `until` (`align`, `cellsnp`, `demux`, `split_bams`, `identify_swaps`), `method` (`solo`, `vireo`, `both`, null), `add_to_previous` (`solo`, `vireo`, null), `hto: multiplexing` (`single`, `multiple`, null). | |
| Warning: Picard runs on STARsolo runs only; `picard_metrics` is ignored for cellranger and cellranger-arc. | Set it to null to silence the warning. |

## Not caught: know these

| Caveat | Why it isn't caught |
|---|---|
| `until` always includes the alignment: the indexed bam, or Picard's metrics when `picard_metrics` is set. Some old `last_step` values (e.g. `STARsolo_demux` without Picard, `*_cellsnp`, `cellranger_demux`, `multiome_gt_demux`) didn't, so their runs now also index the bam (or run Picard). A run that only rebuilds final h5ads from earlier runs (swap correction) needs the alignment outputs to exist; on a project that was aligned they do, and nothing re-runs. | It's how `until` works (`rules/produce_targets.smk`: `UNTIL_STAGES`). |
| `modules` and `references` accept any key, so a typo there (e.g. `modules: STARsolo: bams_dri`) isn't reported; a missing key the rules need fails when the rule is built. | Their keys vary by project (genomes, optional suffixes). |
| `modules: STARsolo: fold_struct` (`Sample_{pool}-cDNA/{pool}-cDNA`) is also the fastq folder of the HTO (kite, with `-cDNA` → `-HTO`) and cDNA fastqs, and the folder of h5ads without demultiplexing: changing it moves those too. | To be split into one `fold_struct` per module (TODO). |
| `demux_run_type` became `genotype: enabled` (`wo_gt` → false, anything else → true). Its other values (e.g. `w_gt_2`) only named the run; `migrate_config.py` says so, and `run_tag` names the run's sub-folder. | The old value isn't in the new config. |
| The unused `fold_struct_gt_demux2` is dropped by `migrate_config.py`. | |
| The schema's `modules: STARsolo: gene_barcodes` default (`_Solo.out/Gene/filtered/barcodes.tsv.gz`) now applies; a duplicate entry without a default overrode it before. It only matters for a config without the key. | |
| The files the config points to (fastqs, `vcf_info` and its columns, VCFs, the wet lab file) are read when the rules are built or run, not by the check. | They are data, not config. |
