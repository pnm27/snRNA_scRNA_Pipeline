# Setting up the Pipeline

## Installing dependencies

### Packages installed through conda

All the packages installed through anaconda3/2018.12 for Python 3.9.5 are described (add_link and file)

### Packages installed through pip

All the python packages installed for python version 3.9.5 ( and with R version 4.1.0) are described (add_link and file)

## Overview of the pipeline

Directed Acyclic Graph of the whole pipeline:
![DAG](images/Whole_pipeline.png)

## Setting up profiles

The info for setting up profiles for different workload managers is mentioned [here](https://github.com/Snakemake-Profiles)

## Executing Pipeline

This pipeline can be executed by executing (in case of any workflow manager, submitting) the script called `run_snakemake.sh`

```sh
sh run_snakemake.sh
```

## Configuration for pipeline

This pipeline depends on a yaml config file (new_config.yaml), which has all relevant options for each rule present in this pipeline. This not only streamlines the process of maintaining or modifying program specific parameters but also makes repetitive usage (can be within the project or using it over mutliple projects) of this pipeline over time very efficient.