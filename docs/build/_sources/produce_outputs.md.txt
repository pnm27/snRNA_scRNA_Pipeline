# Advanced Features

## Configuration File

### Overview

This pipeline depends on a yaml config file (new_config.yaml), which has all relevant options for each rule present in this pipeline. Furthermore, this file has been sectioned in a way that encapsulates rule-specific options/parameters within a sub-workflow module through comments. Each of the rules, however, do have its own required and optional parameters. Thus, necessitating a further demarcation using comments. 

### Common (project-specific) parameters

The following pictures showcase parameters that are project-specific.
:::{figure-md}
![fig1](../../images/config_1.jpg)

new_config.yaml (Part 1)
:::

#### Module selector

**last_step**:
    This is the key which needs to be fed one of the {ref}`pre-selected modules <modules:Selectable Modules>`


