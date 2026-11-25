# Helper scripts

The scripts in `helper_py_scripts/` and `helper_sh_scripts/`, with an example
of each. Some you run yourself, before or after a pipeline run; the others are
run by the pipeline's rules, and their examples are the rules' own calls (to
rerun or debug one step by hand). Every Python script lists all its options
with `--help`; the [API reference](modules) has their functions.

Run the Python scripts from the repository's root, in an environment with
their packages (most need pandas; the h5ad ones need anndata or scanpy, as in
`envs/basic_sctools.yaml`).

## Preparing inputs

### create_wet_lab_info.py: the wet lab info from the lab's spreadsheets

Builds the wet lab info table (`experiment: inputs: wet_lab_info`: the HTOs
and donors of each pool) from the lab's xlsx files. A tracker file lists the
spreadsheets already processed.

```bash
python helper_py_scripts/create_wet_lab_info.py sheets/batch1.xlsx sheets/batch2.xlsx \
    -o projdir/wet_lab_info.tsv -l files_tracker.txt -p PROJ
```

### create_fastq_samp_links_general.bash: links to the fastqs of each sample

Links the fastqs (I1, I2, R1, R2, ...) of each sample into one folder, from a
file listing each sample's R1 fastq (absolute paths). `-n` reads only the first
5 lines.

```bash
find /path/to/fastqs -name "*R1*fastq.gz" | sort > all_r1_fastqs
bash helper_sh_scripts/create_fastq_samp_links_general.bash -i all_r1_fastqs -d projdir/fastq_files/
```

### migrate_config.py: a config in the new layout

Converts a config from the old (flat) layout to the new one, or adds the keys
added since to a config already in the new layout. It lists the changes; the
comments of the example config come along.

```bash
python helper_py_scripts/migrate_config.py old_config.yaml -o new_config.yaml
```

A config already in the new layout, updated in place (back it up first, e.g.
as `backups/new_config.before_<change>.yaml`, which git ignores):

```bash
python helper_py_scripts/migrate_config.py new_config.yaml --in-place
```

### make_wasp_samples.py: a WASP samples table from split bams

Lists per-donor bams as the WASP samples table (`modules: wasp: samples`, for a
`start: wasp` run; see [the WASP samples table](optional_inputs)). Each modality
gets a path pattern with `{pool}` and `{donor}`; bams without a `.bai` are
skipped, and donors with only one modality's bam are reported. `--absolute`
writes absolute paths.

```bash
python helper_py_scripts/make_wasp_samples.py \
    --cDNA 'projdir/STARsolo_split_filt_bams/final/Sample_{pool}-cDNA/{donor}.bam' \
    -o wasp_samples.tsv
```

Multiome, with the ATAC bams split by this pipeline (under `ATAC/`) and the
cDNA bams from elsewhere:

```bash
python helper_py_scripts/make_wasp_samples.py \
    --cDNA '/path/to/rna_split/{pool}/{donor}.bam' \
    --ATAC 'projdir/STARsolo_split_filt_bams/final/ATAC/Sample_{pool}-cDNA/{donor}.bam' \
    -o wasp_samples.tsv
```

### create_filter_perPool_vcf.py and filter_perPool_vcf.bash: filtered per-pool VCFs

Filters each pool's VCF of a folder, as an LSF array job (one task per
VCF). `create_filter_perPool_vcf.py` writes the job script
(`filter_perPool_vcf.bash`, in the current folder) and submits it (`-q`: the
queue); `-n` writes it and lists the VCFs without submitting.
`helper_sh_scripts/filter_perPool_vcf.bash` is a version to submit by hand
(see its `--help`).

```bash
python helper_py_scripts/create_filter_perPool_vcf.py data/vcfs/ results/filtered_vcfs/ -n
```

## After a run

### donor_bams.py: per-donor bams named by the swap correction file

Links (or merges, when a donor is in several pools) the split bams into bams
named by the swap correction file: `<out_dir>/<name>.bam` (`--no_combine`:
`<out_dir>/<name>/<pool>.bam`). Run it with `--dry_run` first: it prints the
plan (the manifest) and writes nothing.

```bash
python helper_py_scripts/donor_bams.py --swap swap_correction.csv \
    --pool_column pool --donor_column vcfID --name_column individualID \
    --bams 'projdir/STARsolo_split_filt_bams/final/Sample_{pool}-cDNA/{donor}.bam' \
    --out_dir donor_bams/ --dry_run
```

### update_logs.py and run_update_logs.sh: one table of every pool's stats

Compiles STARsolo's, Picard's and the demultiplexing stats of each pool,
including STARsolo's run time, into one table. Set the paths in the "CHANGE
THIS" section of `run_update_logs.sh` (`INPUT_FILE`, `BAM_DIR`, `OUTPUT_FILE`,
`BAM_STRUCT`, ...), then:

```bash
bash helper_sh_scripts/run_update_logs.sh
```

### create_aggregate_matrix.py and run_create_aggregate_matrix.sh: one h5ad of all pools

Joins the final h5ads of a folder into one, with the donors renamed by
conversion tables (one per `--suffix`); the top of the script has a worked
example with two tables. `run_create_aggregate_matrix.sh` is an LSF job
running it.

```bash
python helper_py_scripts/create_aggregate_matrix.py final_count_matrix/ Hs_allchr_MT.txt all_pools.h5ad \
    --conversion_file conversion_w_gt.tsv --pool_col Sample --from_don_col Donor_name \
    --to_don_col donor --suffix w_gt --column_prefix SubID_vs_ -g GRCh38
```

### qtltools_plot_samples_wo_gt_multiprocess.py: plots of the genotype checks

Plots QTLtools mbv's matches of each pool's donors (`<pool>_<donor>.bamstat.txt`)
against the vireo donors, with several workers.

```bash
python helper_py_scripts/qtltools_plot_samples_wo_gt_multiprocess.py projdir/mbv/ \
    --qtltools_suffix .bamstat.txt --vireo_out projdir/vireo/ --threads 8
```

### run_scAPAtrap_script.sh: scAPAtrap on a bam

An LSF job running scAPAtrap (an R script, at the path in the job) on a bam,
with samtools and featureCounts from the env modules. Its arguments are passed
on to the R script: the bam, the output folder (its name), the third argument
as is and the fourth as `-t`.

```bash
bsub < helper_sh_scripts/run_scAPAtrap_script.sh   # after editing its #BSUB lines and paths
```

## Run by the pipeline

### create_inp_cellrangerArcCount.py: the libraries csv of cellranger-arc count

Rule `inp_cellranger_arc_count`: sorts a pool's fastqs into gex and ATAC by
their paths.

```bash
python helper_py_scripts/create_inp_cellrangerArcCount.py --output P1_metadata.csv -s P1 \
    projdir/fastq_files/Sample_P1-cDNA/*.fastq.gz projdir/fastq_files/Sample_P1-ATAC/*.fastq.gz
```

### create_Feat_Barc.py: the feature barcodes of a pool's HTOs

Rule `create_FB` (kite): the HTO names and sequences of a pool, from the wet
lab info table (`-c`: its pool, HTO name and HTO sequence columns;
`modules: kite: columns_to_pick`).

```bash
python helper_py_scripts/create_Feat_Barc.py projdir/wet_lab_info.tsv -o P1_FeatureBarcodes.csv \
    -s P1_cDNA -c unique_sample_ID hashtag ab_barcode
```

### create_h5ad_from_bustools.py: an h5ad of the HTO counts

Kite rules: the HTO counts of kallisto/bustools (matrix, genes, barcodes) as an
h5ad.

```bash
python helper_py_scripts/create_h5ad_from_bustools.py output.mtx output.genes.txt output.barcodes.txt \
    -o P1_bustools.h5ad
```

### create_h5ad_from_calico_solo.py: hashsolo's donors of a pool

Rule `run_calico_solo`: runs hashsolo on the HTO counts and keeps the cells
passing the QC (`-m` max mito %, `-g` min genes, `-c` min cells per gene).

```bash
python helper_py_scripts/create_h5ad_from_calico_solo.py P1_bustools.h5ad \
    projdir/STARsolo_outs/Sample_P1-cDNA/Solo.out/GeneFull/raw/ P1_calico_solo.h5ad gene_info.tsv \
    -m 5 -g 1000 -c 10 --mito_prefix MT-
```

### create_inp_cellSNP.py: the barcodes for cellSNP

Rule `create_inp_cellSNP`: the barcodes of a pool's good cells (from its
count matrix, or an h5ad of an earlier run), the input of cellSNP and vireo.

```bash
python helper_py_scripts/create_inp_cellSNP.py projdir/STARsolo_outs/Sample_P1-cDNA/Solo.out/GeneFull/raw/ \
    -o P1_barcodes.txt -m 5 -g 1000 --min_cells 10 --id2name gene_info.tsv --mito_prefix 'MT-'
```

### get_cellSNP_details.awk and get_vireo_details.awk: per-SNP counts

Tables of each SNP's counts: of cellSNP's cells (homozygous and heterozygous
cells) and of vireo's donors (AD and DP per donor). Not run by the rules any
more; run them on the unzipped VCFs.

```bash
zcat cellSNP.cells.vcf.gz | awk -f helper_sh_scripts/get_cellSNP_details.awk > cellSNP_details.txt
zcat GT_donors.vireo.vcf.gz | awk -f helper_sh_scripts/get_vireo_details.awk > vireo_details.txt
```

### demul_samples.py: a pool's h5ad annotated by one demultiplexing tool

Rules `demux_run_<tool>` (Rule A): the pool's count matrix, filtered, with the
donors of hashsolo (`--calico_solo`) or vireo (`--vireo_out`) and their stats.
The rules build the arguments from the config (`lib/params.py`,
`get_params_demux`); a vireo run looks like:

```bash
python helper_py_scripts/demul_samples.py projdir/STARsolo_outs/Sample_P1-cDNA/Solo.out/GeneFull/raw/ \
    P1_vireo.h5ad gene_info.tsv --wet_lab_file projdir/wet_lab_info.tsv \
    -m 5 -g 1000 --min_cells 10 --mito_prefix 'MT-' \
    --vireo_out projdir/vireo/Sample_P1/donor_ids.tsv --stats P1_vireo_stats.tsv --pool_name P1
```

### demux_final.py: a pool's final h5ad

Rule `demux_final` (Rule B): joins the pool's Rule A h5ads (`--sources`: the
tools, in the same order), with the swap correction file's donors when one is
set.

```bash
python helper_py_scripts/demux_final.py --h5ad P1_solo.h5ad P1_vireo.h5ad \
    --stats P1_solo_stats.tsv P1_vireo_stats.tsv --sources 'solo' 'vireo' --pool_name P1 \
    -o P1_final.h5ad --basic_stats P1_basic_stats.tsv --info_out P1_final_info.tsv \
    --swap_correct swap_correction.csv --swap_columns pool donor_from donor_to
```

### collect_demux_versions.py: the demultiplexing versions of the final h5ads

Rule `collect_demux_versions`: one table of the version each pool's final h5ad
was made from.

```bash
python helper_py_scripts/collect_demux_versions.py --pools P1 P2 \
    --basic_stats P1_basic_stats.tsv P2_basic_stats.tsv -o demux_versions.tsv
```

### normalize_external_labels.py: external donor labels as one table

Checkpoint `normalize_external_labels` (`split_by: input: external`): the
cells' donors, from a pooled h5ad or a table, as a pool, Subj_ID, barcodes
table, with the barcodes and pool names transformed to the bams' (JSON lists
of steps, see `transform_steps.py`).

```bash
python helper_py_scripts/normalize_external_labels.py cellBC_vs_origin.txt external_labels.tsv \
    --donor_column individualID --pool_column pool --barcode_column barcode \
    --barcode_steps '[{"type": "pythonic", "operation": "split", "by": "_"}, {"type": "pythonic", "operation": "select_index", "index": -1}]'
```

### create_inp_splitBam_consolidated.py: the barcodes of each donor of a pool

Rule `create_inp_splitBams`: the hash file of a pool (donor, barcode), from an
h5ad, vireo's `donor_ids.tsv`, or the external labels (`--external_pool`).

```bash
python helper_py_scripts/create_inp_splitBam_consolidated.py projdir/vireo/Sample_P1/donor_ids.tsv \
    P1_bc_hash.txt --split_by vireo
```

### check_external_barcodes.py: do the barcodes match the bam?

Rules `check_external_barcodes` and `check_barcodes_atac`: fails when too few
of the bam's CB tags (from its first reads) are the hash file's barcodes, which
is what a barcode format mismatch looks like. `--min_fraction` (default 0.02;
the pipeline passes `modules: split_bams: split_by: external:
min_barcode_fraction`) is low because reads of empty droplets and unlabelled
cells carry CB tags too: a mismatch gives ~0, a pool with few labelled cells
can be at ~5%. To see a pool's fraction on more reads, without failing:

```bash
samtools view pool.bam | head -n 2000000 | \
    python helper_py_scripts/check_external_barcodes.py --hash_file P1_bc_hash.txt --out /dev/null --min_fraction 0
```

```bash
samtools view pool.bam | head -n 200000 | \
    python helper_py_scripts/check_external_barcodes.py --hash_file P1_bc_hash.txt --out P1_bc_check.txt
```

### create_per_donor_bams.bash: a donor's bam

Rules `split_bams` and `split_bams_atac`: the reads of a donor's barcodes,
optionally without the mito reads (`--mito_file`, `--mito_pref`) or outside a
bed (`--bed`), sorted and indexed.

```bash
bash helper_sh_scripts/create_per_donor_bams.bash --donor D1 --hash_file P1_bc_hash.txt \
    --temp_dir projdir/n_temp_d/final/Sample_P1-cDNA/ --out_dir projdir/STARsolo_split_filt_bams/final/Sample_P1-cDNA/ \
    --bam P1_filt.bam
```

### create_per_donor_bams_gt.bash and create_bams_per_donor.sh: older bam splits

Earlier ways of splitting, not run by the rules any more:
`create_per_donor_bams_gt.bash` takes one donor (`samtools view -D`), and
`create_bams_per_donor.sh` merges per-cell bams into donor bams (in batches of
at most the given number of barcodes). Both take positional arguments.

```bash
bash helper_sh_scripts/create_per_donor_bams_gt.bash D1 P1_bc_hash.txt temp/P1/ split/P1/ P1.bam
bash helper_sh_scripts/create_bams_per_donor.sh P1_bc_hash.txt split/P1/ per_cell_bams/P1/ 500
```

### wasp_check_inputs.py: checks of the WASP inputs

The WASP stage's checks before the long jobs. `samples` (rule
`wasp_check_inputs`) checks the samples table against the genotypes;
`contigs` checks that a bam's chromosomes are named as the genotypes';
`set-samples` prints a genotype set's sample IDs (not run by the pipeline: e.g.
to subset your VCFs by hand with `bcftools view -S`).

```bash
python helper_py_scripts/wasp_check_inputs.py samples --table wasp_samples.tsv --set main \
    --names wasp_names.tsv --need_bam -o wasp_inputs_checked.txt
python helper_py_scripts/wasp_check_inputs.py set-samples --table wasp_samples.tsv --set main --sets main --need_bam
python helper_py_scripts/wasp_check_inputs.py contigs --contigs bam_contigs.txt --snp_tab geno/snp_tab.h5
```

### wasp_donor.sh: a donor bam's WASP steps

The steps before (`prepare`, rules `wasp_prepare_rna` / `wasp_prepare_atac`)
and after (`finish`, rules `wasp_finish_rna` / `wasp_finish_atac`) the
remapping, which is a job of its own in between (`wasp_remap_rna`: STAR,
`wasp_remap_atac`: bwa-mem). `prepare` (re)makes the work dir and leaves
`<work>/intersect/to.remap.bam`, `keep.bam` and the fastqs to remap
(`remap.fq1.gz` `remap.fq2.gz`, or `remap.fq.gz`), fixed names whatever the
bam's; `finish` reads them and the remapped bam, and removes the work dir.
Needs the WASP env module (`$WASP_HOME`).

```bash
bash helper_sh_scripts/wasp_donor.sh prepare --work work/P1_D1 --bam D1.bam --samples "G1" \
    --snp_index geno/snp_index.h5 --snp_tab geno/snp_tab.h5 --haps geno/haps.h5 --log D1.intersect.log
bash helper_sh_scripts/wasp_donor.sh finish --work work/P1_D1 --bam D1.bam \
    --remapped work/P1_D1/remapped.bam --out D1.after_wasp.bam \
    --threads 4 --filter_log D1.filter.log --intersect_log D1.intersect.log \
    --flagstat_before D1.before.flagstat --flagstat_after D1.after.flagstat \
    --pool P1 --donor D1 --final_donor D1 --sample G1 --modality cDNA --stats D1.stats.tsv
```

### wasp_run.py: WASP's scripts with a newer numpy

Run by `wasp_donor.sh`: runs one of WASP's python scripts after restoring
the numpy names it uses that newer numpy removed (else: `module 'numpy' has
no attribute ...`): `np.int` (removed in numpy 1.24; `snptable.py`) and
`np.in1d` (removed in numpy 2.4; `find_intersecting_snps.py`). Only names
missing from the numpy in use are set, and numpy's warnings about them are
kept out of the logs. The rest of WASP's scripts runs as it is with Python
3.12 and numpy 2.x (same results as with numpy 1.26). The script's own
arguments follow it.

```bash
python helper_py_scripts/wasp_run.py $WASP_HOME/mapping/find_intersecting_snps.py --is_sorted \
    --snp_index geno/snp_index.h5 --snp_tab geno/snp_tab.h5 --haplotype geno/haps.h5 \
    --samples samples.txt --output_dir work/P1_D1 D1.bam
```

### wasp_allele_counts.py: ref and alt alleles before and after WASP

Run by `wasp_donor.sh finish` with the WASP module's python (it needs pysam
and PyTables, as WASP's own scripts do). Counts the
bases of a donor's reads over its heterozygous SNPs (WASP's own
`snptable.py` and genotype h5 files, as `find_intersecting_snps.py` counts
them) in each bam given: `to.remap.bam` (the reads WASP tested) and
`remap.keep.bam` (those it kept, which are all the reads over SNPs in the
bam after WASP). `--skip_duplicates` leaves out reads marked 0x400 (ATAC).
One row per bam: ref, alt, other, ref_fraction, snps.

```bash
python helper_py_scripts/wasp_allele_counts.py --snp_index geno/snp_index.h5 \
    --snp_tab geno/snp_tab.h5 --haps geno/haps.h5 --samples samples.txt \
    --bam before=to.remap.bam --bam after=remap.keep.bam -o D1.allele_counts.tsv
```

### wasp_stats.py: WASP's statistics

`parse` (run by `wasp_donor.sh finish`) makes a donor bam's stats row from its
logs and flagstats; `summarize` (rule `wasp_stats`, a run of its own) joins a modality's rows
into a table, a summary and plots. `--collect <folder>` adds every donor's
stats file in a modality's folder (`<pool>/stats/*.wasp_stats.tsv`), so the
donors of all batches run into it are joined.

With `--allele_counts` (the table of `wasp_allele_counts.py`), the row has
the reference allele fraction of the reads WASP tested, before
(`summary_ref_allele_fraction_tested`) and after its test
(`summary_ref_allele_fraction_final`, `_change`), and the % of ref and of alt
bases it removed (`summary_ref_removed_pct`, `summary_alt_removed_pct`).
`summary_ref_allele_fraction` is of every read over a SNP before WASP (its
log). The plots: histograms with their kernel density of the `summary_*`
columns that tell donors apart (`PLOT_COLUMNS`: `input`,
`sent_for_remapping_pct`, `failed_remap_pct`, `removed_pct`,
`duplicates_removed_pct`, `ref_allele_fraction`; the rest are counts of
these or their complements, and stay in the table and summary), then a page
of the fraction before vs after, per donor bam. `--log FILE` appends the
time, every stats file joined (with its modification time) and the files
written (rule `wasp_stats` gives each run its own,
`<wasp_dir>logs/wasp_stats.<YYYYmmdd_HHMMSS>.log`).

```bash
python helper_py_scripts/wasp_stats.py parse --pool P1 --donor D1 --final_donor D1 --modality cDNA \
    --sample G1 --intersect_log D1.intersect.log --filter_log D1.filter.log \
    --flagstat_before D1.before.flagstat --flagstat_after D1.after.flagstat -o D1.stats.tsv \
    --allele_counts D1.allele_counts.tsv
python helper_py_scripts/wasp_stats.py summarize P1/stats/*.stats.tsv \
    --stats wasp_stats.tsv --summary wasp_summary.tsv --plots wasp_plots.pdf
python helper_py_scripts/wasp_stats.py summarize --collect wasp_bams/cDNA -o wasp_bams/cDNA \
    --log wasp_bams/logs/wasp_stats.$(date +%Y%m%d_%H%M%S).log
```

## Python modules

Not run directly; imported by the scripts above.

### demultiplex_helper_funcs.py: shared demultiplexing functions

Reading tables, the demultiplexing by hashsolo and vireo, their stats, and the
h5ad annotations, for `demul_samples.py` and `demux_final.py`.

```python
from demultiplex_helper_funcs import auto_read
wet_lab = auto_read("projdir/wet_lab_info.tsv")   # csv, or tsv/txt, by extension
```

### transform_steps.py: the 'pythonic' and 'regex' steps

The transformation steps of the config's recipes (e.g. `split_by: external:
barcode_transformations`), defined in
`schema/demul_samples_annotate_h5ad.schema.json`.

```python
from transform_steps import apply_steps
apply_steps('P12_L1_AAACCTGAGC-1', [
    {'type': 'pythonic', 'operation': 'split', 'by': '_'},
    {'type': 'pythonic', 'operation': 'select_index', 'index': -1},
])  # 'AAACCTGAGC-1'
```
