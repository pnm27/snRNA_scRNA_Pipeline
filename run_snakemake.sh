#!/bin/bash

#BSUB -J Gen_Snakemake
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=1000]
#BSUB -W 48:00
#BSUB -oo /sc/arion/projects/psychAD/pnm/Snakemake.stdout
#BSUB -eo /sc/arion/projects/psychAD/pnm/Snakemake.stderr
#BSUB -L /bin/bash


# For snakemake < v8
# source /sc/arion/work/prashf01/conda/envs/snakemake/etc/profile.d/conda.sh
# conda activate snakemake


# snakemake --profile prachu_lsf --restart-times 2 --latency-wait 10 #--batch all=1/3

# For snakemake > v8
source /sc/arion/work/prashf01/conda/envs/new_snakemake/etc/profile.d/conda.sh
conda activate new_snakemake
export SNAKEMAKE_LSF_MEMFMT="perjob"
snakemake --profile prachu_lsf_smk8 --workflow-profile workflow_profile \
 --rerun-triggers input software-env mtime # code params
