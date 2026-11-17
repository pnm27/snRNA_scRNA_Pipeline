# Overview of the Pipeline
To set up the pipeline for projects, the following setups need to be performed:

1. Setting up wildcards for the project (see {ref}`how <sub_workflows:Wildcard Processing>`).
2. According to the experiment select a module (module options).
3. As required and wherever possible, setup the folder structures for different programs.
4. Check and modify parameters of the programs to be run (in new_config.yaml file).
   1. (Currently) make sure all required python and R packages are present.
5. Change global variables as required (in Snakefile - check).

# Wildcard Processing
For the purpose of creating wildcards a list of samples to be processed is provided to the pipeline. There are 3 ways to achieve this:

- list of samples/pools (as a folder structure)
- yaml file containing the list of samples/pools
- Directory containing input files

## List of samples



``` {currentmodule} helper_py_scripts/

```

This pipeline has many combinations of the aforementioned programs as a built-in set that can be executed using specific keywords.

```{eval-rst}
.. _selectable-modules:
```

# Selectable Modules
What a run does is set by `workflow: aligner`, `workflow: start` and `workflow: until` in the config:

- `aligner`: `STARsolo`, `cellranger` (cellranger count) or `cellranger_arc` (cellranger-arc count, multiome)
- `start`: the first stage to run: `align` (from the fastqs) or `split_bams` (split bams that are already there by donor labels from outside the pipeline; see [Splitting bams by external donor labels](#splitting-bams-by-external-donor-labels))
- `until`: the last stage to run. The alignment is always a target (unless `start: split_bams`): the indexed bam, or Picard's metrics when `workflow: picard_metrics` is set (STARsolo only).

| `until` | Stages run |
|---|---|
| `align` | alignment (or Picard's metrics) |
| `cellsnp` | alignment, cellSNP's pileup |
| `demux` | alignment, demultiplexing (`experiment: demultiplexing: method`) |
| `split_bams` | alignment, demultiplexing, per-donor bams |
| `identify_swaps` | alignment, demultiplexing, per-donor bams, QTLtools-mbv |

Multi-VCF runs (one vireo run per vcf column of `vcf_info`) set `experiment: demultiplexing: genotype: multi_vcf: true`.

`helper_py_scripts/migrate_config.py` converts configs with the earlier `last_step` (e.g. `STARsolo_gt_demux_identify_swaps` → `aligner: STARsolo`, `until: identify_swaps`). The module names in the table below are those `last_step` values.

## Module description

```{list-table} Modules_info
:header-rows: 1
:name: module-info-table

* - Module Name
  - Module Info
  - Sub Worflows Involved
* - all
  - module_info (more desc in its own file)
  - sub_wkfl
* - all_multi_vcf
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_kb_solo
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gt_demux
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_split_bams
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_split_bams_gt_demux
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_split_bams_gt_demux_multi_vcf
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gt_demux_multi_vcf
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_cellsnp
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_gt_demux_identify_swaps
  - module_info (more desc in its own file)
  - sub_wkfl
* - starsolo_resolve_swaps_gt_demux
  - module_info (more desc in its own file)
  - sub_wkfl
```

# Sub-Snakemake workflows
This pipeline divides each module into its self-contained individual workflows. These are:

```{list-table} Sub_Snakemake_workflows_table
:header-rows: 1
:name: sub_smk-info-table

* - Name of Workflow
  - Description
* - resources.snkmk 
  - It contains memory (in MB per thread) and time requirements (in minutes) for each rule.
* - calico_solo_demux.snkmk 
  - It contains hashsolo rule.
* - split_bams.snkmk [^consol]
  - It contains rules needed to split pooled bams into individual bams dependent on output
  produced by either hashsolo or vireoSNP using custom scripts.
* - input_processing.snkmk 
  - It contains rules that collects values for all the wildcards.
* - STARsolo.snkmk 
  - It contains rules for STARsolo.
* - produce_targets.snkmk 
  - It contains the **rule all** and the needed functions.
* - snv_aware_align.snkmk [^not-yet]
  - **This might be removed soon**
* - kite.snkmk  
  - It contains rules for the [kite](https://github.com/pachterlab/kite#kite-kallisto-indexing-and-tag-extraction) workflow.
* - picard_metrics.snkmk 
  - It contains rules for all PICARD metrics (GCBiasMetrics and RNAseqMetrics).
* - pheno_demux3.snkmk 
  - It contains rules for the cellSNP-vireoSNP pipeline.
* - split_bams_gt.snkmk [^consol]
  - It contains rules needed to split pooled bams into individual bams dependent on output produced by vireo.
* - demultiplex.smk
  - It contains the rules that create an h5ad per demultiplexing tool (`demux_run_<tool>`) and the final h5ad of each pool (`demux_final`). See [Demultiplexing outputs](#demultiplexing-outputs).
* - identify_swaps.snkmk 
  - It contains rules for identifying swaps using QTLtools-mbv.
```

[^consol]: Consolidating into one

[^not-yet]: Not yet implemented

# Demultiplexing outputs

Each demultiplexing run writes to its own sub-folder (`experiment: demultiplexing: run_tag` in the config, e.g. `wo_gt` or `w_gt_2`) of the calico_solo, cellSNP, vireoSNP, final count matrix and demultiplex info dirs. With `run_tag: null` the dirs are used as they are.

```{list-table} Demultiplexing rules
:header-rows: 1
:name: demux-rules-table

* - Rule
  - Outputs
* - `demux_run_solo`, `demux_run_vireo` (`demux_run_none` without demultiplexing)
  - An h5ad with the classification of one tool (QC stats in `uns['qc_stats']`) and the stats of that tool (`_hashsolo_stats.tsv` / `_vireo_stats.tsv`), in the tool's dirs.
* - `demux_final`
  - The final h5ad of each pool (`modules: demultiplex: demux_count_matrix_dir`, `modules: demultiplex: fold_struct`), with the donor of each cell in `Donor_final`; its basic stats (`_basic_stats.tsv`: QC and the demultiplex versions used) and its demultiplex info (basic and tool stats together, read by `update_logs.py`).
* - `collect_demux_versions`
  - `final_demux_versions.tsv` in `modules: demultiplex: demultiplex_info_dir`: one row per pool with the demultiplex versions used, whether the swap correction was applied, and the number of donors and of cells passing QC.
```

`demux_final` combines, for each pool:

1. the runs in its row of the swap correction file (`modules: demultiplex: swap_correction`), given as `<tool>:<run_tag>` (tool: `solo`, `vireo` or `both`) or as a name in `versions`, whose `donor_from` → `donor_to` columns rename the donors (an empty `donor_to`, or a `status_column` value other than True, marks the donor's cells `Removed`), else
2. the tools of this run and, with `add_to_previous`, the other tool from the run `modules: demultiplex: previous_run_tag`.

To split the bams by the final h5ad (e.g. after the swap correction), set `modules: split_bams: split_by: input: final`.

## Splitting existing bams by the versions of the swap correction file

When the bams and the demultiplexing runs already exist, a run can start at the split (`workflow: start: split_bams`) and take everything from the swap correction file, e.g.

```text
Sample   Donor_name  Solution  perc_het_consistency  ...  final_status  Version_used
Sample1  donor2      donA      0.953693              ...  True          w_gt
Sample1  donor3      donB      0.955703              ...  True          w_gt
Sample2  donor4      donC      0.960862              ...  True          w_gt
```

```yaml
workflow:
  start: split_bams
  until: split_bams
modules:
  demultiplex:
    swap_correction:
      file: swap_corrections.tsv
      pool_column: Sample
      version_column: Version_used
      donor_from: Donor_name
      donor_to: Solution
      status_column: final_status
      versions: demul_samples_annotate_h5ad_template.json
  split_bams:
    split_by:
      demux: vireo
      input: final
      external:
        pools: bams.tsv   # optional: bams outside the aligner's folders
```

- The pools are those of the file with a donor whose `status_column` is True (`True`, `yes` or `1`); the others aren't split (a warning). Donors not True are dropped (`Removed`), as an empty `donor_to` is.
- `versions` says what the names in `version_column` are: a JSON file whose `demultiplex_paths` gives each name's `demux` (h5ad) dir and, optionally, `info` (stats) dir, as in the h5ad annotation template (so the template and the run share them); or a map in the config, of each name to `<tool>:<run_tag>` (a run of this pipeline), or to `{demux, info, tool}`. Without it, the column holds `<tool>:<run_tag>`.
- In a version's dirs the file names follow the tool's `fold_struct` and `final_count_matrix_h5ad` (vireo: `<demux>Sample_<pool>-cDNA_STARsolo_vS_out.h5ad` and `<info>Sample_<pool>-cDNA_vireo_stats.tsv` with the example's values); the tool is `split_by: demux` unless the map gives one.
- Only `demux_final` (the final h5ads, in `modules: demultiplex: demux_count_matrix_dir`) and the split run. The bams are in the aligner's folders, or in the pools table (`split_by: external: pools`), as with external labels; pools of the swap correction file missing from the table aren't split (a warning).
- Keys left out take the defaults of `schema/split_bams_defaults.yaml`, as below.

Split bams go to a sub-folder of `inp_split_bams_dir`, `split_bams_dir`, `per_donor_split_log_dir` and `new_temp_dir`: `gt_check/` for genotype checks (`workflow: gt_check: yes`; the whole bam, or the one subset to `subset_chr`) and `final/` otherwise (the bam filtered by cell barcode).

## Splitting bams by external donor labels

Bams can be split by donor labels made outside the pipeline, e.g. a pooled h5ad annotated after the analysis of many pools. The file has, for each cell, its donor and its pool of origin:

- a pooled h5ad: obs columns (the barcodes are the obs names, unless `barcode_column` is set), or
- a table (`.csv`, or tab-separated `.tsv` / `.txt`; optionally `.gz`), one row per cell:

```text
barcode                     individualID   pool_of_origin
NPSAD-P12_L1_AAACCTGAGC-1   R1234567       P12
```

To split the bams the pipeline aligned earlier (in the aligner's `bams_dir`, named by its `fold_struct`), without aligning or demultiplexing again:

```yaml
workflow:
  aligner: STARsolo
  start: split_bams
  until: split_bams     # or identify_swaps
  gt_check: no          # yes: bams for genotype checks

modules:
  split_bams:
    split_by:
      input: external
      external:
        file: projdir/pooled_annotated.h5ad
        donor_column: individualID
        pool_column: pool_of_origin
        barcode_column: null
        # From the file's unique barcodes to the bam's CB tags
        barcode_transformations:
          - {type: pythonic, operation: split, by: "_"}
          - {type: pythonic, operation: select_index, index: -1}
        # From the file's pool names to the bams' {pool}, e.g. P12 -> NPSAD-P12
        pool_transformations:
          - {type: regex, pattern: "^", replacement: "NPSAD-"}
```

With `start: align`, the pools are aligned first and then split by the labels (no demultiplexing). The transformation steps are those of the h5ad annotations (`schema/demul_samples_annotate_h5ad.schema.json`): `pythonic` (`lower`, `upper`, `title`, `strip`, `rstrip`, `lstrip`, `split` / `join` with `by`, `select_index` with `index`, `select_range` with `start` / `end`, `replace` with `old` / `new`) and `regex` (`pattern`, `replacement`, optional `flags`). E.g. to add cellranger's `-1`: `{type: regex, pattern: "$", replacement: "-1"}`.

The jobs:

1. `normalize_external_labels` (checkpoint) reads the labels once, applies the steps and drops doublets, negatives and unassigned cells: `inp_split_bams_dir/external_labels.tsv` (pool, Subj_ID, barcodes). Its log lists the cells and donors of each pool and example barcodes. Pools without a bam stop the run (`start: split_bams`).
2. `create_inp_splitBams` writes each pool's donor hash file (`donor_name_converter` and `min_cells_per_donor` apply).
3. `check_external_barcodes` stops the run if the hash file's barcodes don't look like the bam's CB tags, before any bam is filtered or split.
4. `bamfilt_by_CB` / `filt_chr_bams` and `split_bams` as for the other inputs.

### Bams outside the aligner's folders (pools table)

When the bams aren't where the aligner's `bams_dir` and `fold_struct` put them, e.g. `root/dayA/samp1/`, `root/dayB/samp3/`, list them in a pools table (`split_by: external: pools`, with `start: split_bams`). Columns: `pool`, `bam` and, optionally, `bai` (default: `<bam>.bai`); tab-separated, or comma-separated for `.csv`:

```text
pool    bam
samp1   /abs/root/dayA/samp1/samp1_Aligned.sortedByCoord.out.bam
samp3   /abs/root/dayB/samp3/samp3_Aligned.sortedByCoord.out.bam
```

```bash
{ printf 'pool\tbam\n'; for b in root/day*/samp*/*_Aligned.sortedByCoord.out.bam; do printf '%s\t%s\n' "$(basename "$(dirname "$b")")" "$(realpath "$b")"; done; } > bams.tsv
```

```yaml
workflow:
  start: split_bams
  until: split_bams
modules:
  split_bams:
    split_by:
      input: external
      external:
        file: projdir/cellBC_vs_origin.txt
        pool_column: origin
        barcode_column: barcode
        pools: projdir/bams.tsv
```

That can be the whole config: with `start: split_bams`, keys left out take the values of `new_config.example.yaml`, with `schema/split_bams_defaults.yaml` on top (no demultiplexing, Picard metrics or genotype check; `aligner: STARsolo`; `donor_column: individualID`; the barcode steps keep the part after the last `_`). No fastq list is read. Set what differs, e.g. the output dirs of `modules: split_bams` (default `projdir/...` in the run's dir), `barcode_transformations: []` when the labels' barcodes are the CB tags, or `gt_check: yes` with `gt_demux: vcf_info`. Runs from the fastqs (`start: align`) get no defaults.

- The `pool` names are the labels' pools after `pool_transformations`; they can't contain `/`. Pools of the labels missing from the table aren't split, with a warning (e.g. when only some bams were kept, for disk space); the run stops when none of them is in the table. Pools in the table whose bam or bai is missing stop the run (listed); pools only in the table aren't split (noted once).
- The bams are only read: the cell-filtered and chromosome-subset bams (`filt_bam`, `short_bam`) go to `split_bams: temp_dir` (`<temp_dir><pool>/<pool>_filt.bam`) instead of next to the bam.
- `workflow: aligner` still sets the barcode format (cellranger's `-1`); multiome (`cellranger_arc`) isn't supported with a table.

## Per-donor bams named by the swap correction file (stand-alone)

`helper_py_scripts/donor_bams.py` (not a pipeline rule yet, see TODO) takes split per-donor bams and names them by the swap correction file: each row's pool and donor (`--donor_column`, the name of the split bam, e.g. the vcf ID) get a name (`--name_column`, e.g. individualID). A name with one bam is linked to `<out_dir>/<name>.bam` (hard link; a symlink across filesystems, `--copy` to copy); a name with bams of several pools is merged (`samtools merge`) into one, each read with a read group of its pool (`RG:Z:<pool>`, as the same cell barcode in two pools is two cells). Donors with an empty name, or a `--status_column` value other than True / yes / 1, are dropped.

```bash
python3 helper_py_scripts/donor_bams.py --swap final_swap_correction.csv \
    --pool_column pool --donor_column vcf_id --name_column individualID \
    --status_column final_status \
    --bams 'projdir/STARsolo_split_filt_bams/final/Sample_{pool}-cDNA/{donor}.bam' \
    --out_dir projdir/donor_bams/ --threads 8 --dry_run
```

- `--keep_apart tissue`: merge only the bams that also share these columns (the same donor from different tissues or conditions), as `<name>_<tissue>.bam`
- `--no_combine`: never merge; each bam is linked to `<out_dir>/<name>/<pool>.bam`
- `--dry_run` prints the plan; otherwise it is written to `<out_dir>/manifest.tsv` (output, action: link / copy / merge / dropped, source bams, pools, donors). Missing bams are all listed before anything is written; bams with different references (`@SQ`) aren't merged; outputs that exist are kept unless `--force`

## WASP filtering of per-donor bams

For sc/snRNA-seq and the RNA (GEX) part of multiome (STARsolo, cellranger or cellranger-arc bams; ATAC isn't covered yet, see TODO). Each per-donor bam is filtered with [WASP](https://github.com/bmvdgeijn/WASP)'s mapping pipeline: reads over the donor's heterozygous SNPs are remapped with the alleles swapped and kept only if they map back to the same place (no allele-specific mapping bias), along with the reads over no SNP. Duplicates are kept (sc/snRNA-seq). This is different from STARsolo's own WASP mode (`modules: STARsolo: wasp`, the `vW` tag at alignment).

Two ways to run it:

- `workflow: until: wasp`: after `split_bams`, on this run's split bams (`gt_check: no`: the final, cell-filtered bams). It follows any `split_bams` input, including external labels and `start: split_bams`.
- `workflow: start: wasp` (with `until: wasp`): on per-donor bams that are already there, listed in the samples table. No alignment or demultiplexing; keys left out take the defaults of `new_config.example.yaml` with `schema/wasp_defaults.yaml` on top.

### The samples table

`modules: wasp: samples`, tab-separated (or comma-separated for `.csv`), one row per donor bam:

| column | |
|---|---|
| `pool` | pool (no `/`) |
| `donor` | donor (no `/`); names the output bam. After `split_bams`, the split bam's donor (e.g. its vcf ID) |
| `sample` | the donor's genotype sample ID(s) in its genotype set, separated by `;` (or `,` in a tab-separated table), no whitespace in an ID: written one per line to the `--samples` file of `find_intersecting_snps.py`. Optional: empty, or no column, takes the ID from `modules: wasp: sample_names`, else the donor name |
| `bam`, `bai` | `start: wasp` only: the per-donor bam (`bai` defaults to `<bam>.bai`) |
| `genotype_set` | which `modules: wasp: genotypes` set, with more than one |

```text
pool    donor   sample      bam
P1      D1      GENO_0001   /data/split/Sample_P1-cDNA/D1.bam
P1      D2      GENO_0002   /data/split/Sample_P1-cDNA/D2.bam
```

After `split_bams`, split donors not in the table aren't filtered (a warning); none at all stops the run. With `start: wasp`, a missing bam or bai stops the run, all listed.

When pools or donors aren't named as in the genotypes (e.g. one donor per pool, named after the pool), `modules: wasp: sample_names` maps them to genotype IDs instead of filling `sample` row by row: columns `pool`, `sample` and optionally `donor` (empty: every donor of the pool). See [Optional input files](optional_inputs.md#wasp-sample-names).

### Genotypes

Named sets (e.g. one for WGS, one for SNP arrays), each:

- `h5_dir`: a folder with WASP's `snp_index.h5`, `snp_tab.h5` and `haps.h5` (made once with WASP's `snp2h5`), or
- `vcf` (an indexed VCF/BCF, or a list of them, e.g. one per chromosome) and `chrom_sizes` (chromosome and length, `snp2h5 --chrom`): rule `wasp_snp2h5` makes the h5 files once (`<wasp_dir>genotypes/<set>/`), from the set's samples in the table and the sites variable among them (`bcftools view -S ... -q 0.000001:minor`), one chromosome at a time.

### Configuration

```yaml
references:
  tools:
    wasp: wasp/<version>     # env module that sets WASP_HOME (its scripts are under it)
modules:
  wasp:
    samples: wasp_samples.tsv
    genotypes:
      main:
        h5_dir: /path/to/snp2h5_output/
      # chip: {vcf: /path/to/genotypes.bcf, chrom_sizes: /path/to/chromInfo.hg38.txt}
    remap:
      genome: null           # STAR index; null: the aligner's reference (start: wasp: required)
      gtf: null              # null: STARsolo's genome_pick gtf when the genome is the aligner's
      extra_params: null     # more STAR options
    paired_end: false        # true for paired-end GEX reads; 10x 3': false
    wasp_dir: projdir/wasp_bams/
    work_dir: projdir/wasp_tmp/   # large temporary files, removed per job: a place with space
```

The remap genome must be the one the bams were aligned to: by default STARsolo's `genome_pick` (with its gtf), or the `star` folder of the cellranger / cellranger-arc reference (see `config_caveats.md`).

### Jobs

1. `wasp_snp2h5` (sets given as VCF/BCF only).
2. `wasp_check_inputs`, once: every `sample` of the table is in its set's `haps.h5`, so a wrong ID stops the run before the long jobs.
3. `wasp_remap`, per donor bam (as the lab's tested scripts, e.g. `genesis_wasp.R`): a check that the bam and the genotypes name chromosomes the same way (`chr1` vs `1`); `find_intersecting_snps.py`; STAR remapping (two-pass, `NH HI AS nM NM MD`); `filter_remapped_reads.py`; `samtools merge`, `sort`, `index`; `samtools flagstat` of the bam before and after; the donor's stats. Each job holds one `wasp_slots`: `workflow_profile/config.yaml` sets `resources: wasp_slots=100`, so at most 100 run at a time (`--resources wasp_slots=N` on the command line to change it). 6 cores, 18 GB, a day by default (`execution: resources: wasp_remap`).
4. `wasp_stats`: the run's table, summary and plots.

Outputs, under `wasp_dir`:

```text
<pool>/<donor>.after_wasp.bam(.bai)
<pool>/stats/<donor>.wasp_stats.tsv, .before_wasp.flagstat, .after_wasp.flagstat
<pool>/logs/<donor>.find_intersecting_snps.log, .star_remap.log, .filter_remapped_reads.log
wasp_stats.tsv     one row per donor bam
wasp_summary.tsv   mean, median, sd, min, max of each summary_* column
wasp_plots.pdf     histograms of the summary_* columns
```

The stats (`helper_py_scripts/wasp_stats.py`, as in `psychad_multiome_wasp.R`) are the reasons reads were set aside by `find_intersecting_snps.py` (`intersect_discard_*`), its kept and to-remap reads and SNP allele matches, `filter_remapped_reads.py`'s counts (`filter_*`), and the flagstat counts. The `summary_*` columns: `input` (mapped reads before), `no_snp_kept`, `sent_for_remapping`, `failed_remap` (bad + discarded), `final` (mapped reads after), `removed`, each also as `_pct` of `input`, and `ref_allele_fraction` of the reads over SNPs. Counts are reads (pairs with `paired_end: true`).

### Usage

WASP only, from a samples table with bams:

```yaml
workflow:
  start: wasp
  until: wasp
references:
  tools:
    wasp: wasp/<version>
modules:
  wasp:
    samples: wasp_samples.tsv
    genotypes:
      main:
        h5_dir: /path/to/snp2h5_output/
    remap:
      genome: /path/to/STAR_index/
      gtf: /path/to/genes.gtf
    work_dir: /path/with/space/wasp_tmp/
```

```bash
snakemake --configfile wasp_config.yaml --profile <cluster profile> --workflow-profile workflow_profile --use-envmodules --use-conda
```

After a run that split bams, the same config with `until: wasp` and a samples table without the bam column (`--rerun-triggers mtime` keeps the split bams from being remade when only the config changed).

Several datasets in one table, summary and plots:

```bash
python3 helper_py_scripts/wasp_stats.py summarize A4/wasp_bams/wasp_stats.tsv A13/wasp_bams/wasp_stats.tsv \
    --labels A4 A13 -o wasp_all/
```

```{eval-rst}
.. autosummary::
  :nosignatures:

   demultiplex_helper_funcs._ret_htos_calico_solo  
```
