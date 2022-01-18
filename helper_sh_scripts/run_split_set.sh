#!/bin/bash

#BSUB -J split_h5ad_files
#BSUB -P acc_CommonMind
#BSUB -q express
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=5000]
#BSUB -W 03:00
#BSUB -oo /sc/arion/projects/psychAD/pnm/split_h5ad_files.stdout
#BSUB -eo /sc/arion/projects/psychAD/pnm/split_h5ad_files.stderr
#BSUB -L /bin/bash


ml anaconda3
ml R/4.1.0

source /sc/arion/work/prashf01/conda/envs/snakemake/etc/profile.d/conda.sh
conda activate snakemake

python3 /sc/arion/projects/psychAD/pnm/helper_py_scripts/split_seth5ad_to_samph5ad.py
