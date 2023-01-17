% scRNAseq Pipeline documentation master file, created by
% sphinx-quickstart on Thu Dec 29 16:11:21 2022.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

```{warning}
This documentation is incomplete and is under heavy development!
```

```{eval-rst}
.. admonition::
   TODO:

   * Add rule-specific `conda envs <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management>`_ and verify them.
      * Conda envs can't be used with **run**. Hence, appropriately modify rule using **run** directive to facilitate usage of conda env or make a note of it in the documentation.
   * Move documentation to configargparse.
   * Write down schemas.
   * More documentation.
   * Combine rules split_bams and split_bams_gt.
   * Add new Picard metrics.
   * Include in new_config.yaml an option to select wasp mode in the rule STARsolo
   * Search Ranking of readthedocs (using config file for this too).
   * Might think to incorporate git submodules for repos on git that I use.

```

```{include} ../../README.md
:relative-images:
```

```{toctree}
:caption: Introduction
:maxdepth: 2

highlights
set_up
```

```{toctree}
:caption: Sub Workflows
:maxdepth: 2

modules
```

```{toctree}
:caption: Advanced Setups
:maxdepth: 3

produce_outputs
```

```{toctree}
:caption: Reference
:maxdepth: 3

glossary
```

# Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
