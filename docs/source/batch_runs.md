# Running in batches

A long list of pools or donors can be run a part at a time, e.g. to stay within a disk quota, to start while some inputs are still being made, or because some bams were moved.

Every stage takes its pools (or donors) from one list file. A batch is that list with fewer rows. The output dirs stay the same for every batch, so the batches add up to the same files one run of the whole list makes.

| Run | List that sets the batch | Per batch, change | Keep the same |
|---|---|---|---|
| [From the fastqs](#from-the-fastqs-align-cellsnp-demux) (`start: align`) | fastq list (`experiment: inputs: select_fastqs`) | the rows of the list | everything else |
| [Split by external labels](#split-by-external-labels) (`split_by: input: external`) | pools table (`split_by: external: pools`) | the rows of the table | the labels file, the output dirs |
| [Bams moved](#bams-moved-to-another-location) | pools table | the `bam` of the moved pools | the labels file, the output dirs |
| [Split by the swap correction file](#split-by-the-swap-correction-file) (`split_by: input: final`) | pools table, else the swap correction file | the rows of the table | the versions, the output dirs |
| [WASP](#wasp) (`start: wasp`) | samples table (`modules: wasp: samples`) | the rows of the table | `wasp_dir` |

Either edit the list between batches, or keep one list per batch and a config for each that differs only in the path of the list (as in the WASP example, "Samples in batches" under [Usage](sub_workflows.md#usage)).

## From the fastqs (align, cellSNP, demux)

**Change:** the [fastq list](optional_inputs.md#fastq-list). The pools of a run are its lines, so a batch is the list with some lines, or with a `#` in front of the pools to leave for later:

```text
Sample_P1-cDNA
Sample_P2-cDNA
# Sample_P3-cDNA
# Sample_P4-cDNA
```

**Leave as they are:** the pool info file (`vcf_info`), the wet-lab sheet and the swap correction file. They are looked up by pool, so they can hold every pool from the first batch on; rows of pools not in the run aren't read.

**Why the list:** it is the only place the pools of a `start: align` run come from. Outputs of the pools of earlier batches aren't targets of a later one, so they are neither checked nor redone.

With final h5ads (`until: demux` or later), the table of the demultiplexing versions used (`collect_demux_versions`) is one file for the run's pools, so each batch writes it with its own pools. After the last batch, run once with every pool in the list: the pools' files are there, so only that table is left to make.

## Split by external labels

The pools to split are those of the labels (`split_by: external: file`) that have a bam. There are two lists to cut, and the pools table is the one to prefer.

### By the pools table (recommended)

**Change:** the [pools table](optional_inputs.md#pools-table-bamstsv) (`split_by: external: pools`): some of its rows in each batch.

```text
pool    bam
samp1   /abs/root/dayA/samp1/samp1_Aligned.sortedByCoord.out.bam
```

**Leave as they are:** the labels file and `inp_split_bams_dir`, `split_bams_dir` and the other dirs of `modules: split_bams`.

**Why:** pools of the labels that aren't in the table aren't split (a warning lists them), so a shorter table is a batch. The labels are read once (`normalize_external_labels`, which writes `<inp_split_bams_dir>external_labels.tsv`) and that table stays valid for every batch, so no batch redoes it.

- A pool listed again in a later batch isn't redone: its files are there and up to date. So the table can also grow (batch 1's rows, then batch 1's and 2's) instead of being replaced.
- A pools table needs `workflow: start: split_bams`. When the bams are in the aligner's folders (no pools table so far), write one for them: `pool` and the path the aligner gave each bam. Note that with a table the cell-filtered bams go to `split_bams: temp_dir` instead of next to the bam.
- Multiome: the rows also hold `atac_bam`, so a pool's gex and ATAC bams are in the same batch.

### By the labels file

**Change:** `split_by: external: file`, to a file with the cells of the batch's pools only (a table, or an h5ad subset by its pool column).

**Leave as they are:** the pools table, if any (pools only in the table aren't split), and the output dirs.

**Why it's the second choice:** each batch's file is newer than `external_labels.tsv`, so each batch reads the labels again and writes that table with its own pools only. The pools of earlier batches aren't redone, as they aren't in it. But a later run of all the labels, e.g. the whole file to check that nothing was left out, finds `external_labels.tsv` newer than every pool's donor hash file and splits every pool again. Use it when a pools table isn't wanted, e.g. with `start: align`, where the fastq list has to be cut to the same pools too.

`external_labels.tsv` itself is the pipeline's output: change the file of `split_by: external: file`, not that table.

## Bams moved to another location

**Change:** the `bam` column (and `bai`, `atac_bam`, `atac_bai`, if listed) of the moved pools in the pools table. If the run had no pools table, because all bams were in the aligner's folders, add one with every pool still to split and set `split_by: external: pools` and `workflow: start: split_bams`.

**Leave as they are:** the labels file and the output dirs.

**Why:** with a pools table each pool's bam is read where its row says, and nothing else in the config holds bam paths. Pools whose bam or bai isn't at the listed path stop the run before any job, with the list of them, so a path that wasn't updated is caught at the start.

- Pools not yet split are simply read from the new place.
- Pools already split have a new input path. Snakemake redoes jobs whose set of input files changed, so they would be split again; leave them out of the table, or run with `--rerun-triggers mtime` (moving a file keeps its time, so they are then seen as up to date).
- Bams that were removed (e.g. for disk space): leave their rows out. Their pools are named in the warning and the rest is split.

## Split by the swap correction file

`workflow: start: split_bams` with `split_by: input: final`: the pools are those of the [swap correction file](optional_inputs.md#swap-correction-file) with a donor whose `status_column` is True.

**Change:** the pools table (`split_by: external: pools`), as for external labels: its pools that are missing from the table aren't split (a warning). Without a pools table, the swap correction file is the list: give each batch a copy with the rows of its pools.

**Leave as they are:** `versions` and the dirs of the demultiplexing runs, the output dirs.

**Why the pools table first:** the swap correction file is also read for the final h5ads and the donor names, so one file for all batches keeps those the same; the pools table only says which bams are there.

## WASP

**Change:** the samples table (`modules: wasp: samples`): some rows per batch. **Keep:** `wasp_dir`. The stats of all donors there are joined by a run of their own after the batches (`snakemake --configfile <any batch's config> wasp_stats --use-conda`), not by the batches; each such run's log (`<wasp_dir>logs/wasp_stats.<time>.log`) lists the donors' stats files it joined. See "Samples in batches" under [Usage](sub_workflows.md#usage) of the WASP section.

## For every batch

- **One after another in a folder.** Snakemake locks the folder a run is started from, so a second batch started there at the same time stops with "Directory cannot be locked". Batches started from different folders that write to the same output dirs aren't protected from each other: give them pools that don't overlap, and don't let two of them read the labels at once (both would write `external_labels.tsv`).
- **A pool in one batch.** Listing a pool again is harmless for the split (it isn't redone), but a WASP donor listed again is redone.
- **Dry run first.** `snakemake -n` with the batch's config lists the jobs and prints the warnings of pools left out, which should be the pools of the other batches and no others. With external labels that haven't been read yet, the dry run stops at `normalize_external_labels`, as the pools aren't known before it.
- **Same output dirs and `run_tag`.** A batch with other dirs or another `run_tag` is another run, and its files don't add to the earlier batches'.
