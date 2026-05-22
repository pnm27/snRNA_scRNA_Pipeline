# snRNA_scRNA_Pipeline Introduction

## Why the need?

This pipeline intends to not only make complex {term}`preprocessing` workflows easy (e.g. snRNA seq with pooled samples, double HTOs, etc.) but also to facilitate the use of common workflows used for preprocessing by providing *readymade* different combinations of softwares/tools (see {ref}`selectable <selectable-modules>` modules for more options).

It also supports various software/pipeline for scRNA seq pre-processing.

The highlights of the pipeline are:

* Streamlined processes to modify parameters for each program through a single yaml file
* Easily modifiable to accomodate more rules
* Can be used for both individual samples as well as multiplexed pools
* Preserve folder structures (mirroring fastqs' folder structures)
* Organize outputs from each module
* Select multiple pre-set modules that simplifies usage across multiple projects

## Software Used

This pipeline depends on the following packages/programs:

- [Scanpy Manual](https://scanpy.readthedocs.io/en/stable/)
- [Snakemake Manual](https://snakemake.readthedocs.io/en/stable/)
- [STARsolo Manual](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)
- [GC bias metrics (PICARD)](https://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics)
- [RNA seq metrics (PICARD)](https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics)
- [KITE *(kallisto indexing and tag extraction)*](https://github.com/pachterlab/kite)
- [cellSNP Manual](https://cellsnp-lite.readthedocs.io/en/latest/manual.html)
- [vireoSNP Manual](https://vireosnp.readthedocs.io/en/latest/manual.html)
- [hashsolo Info](https://github.com/calico/solo#how-to-demultiplex-cell-hashing-data-using-hashsolo-cli)
- [QTLtools-mbv](https://qtltools.github.io/qtltools/pages/QTLtools-mbv.1.html)
