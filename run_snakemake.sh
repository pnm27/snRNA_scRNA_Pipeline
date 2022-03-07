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


ml anaconda3
ml R/4.1.0

source /sc/arion/work/prashf01/conda/envs/snakemake/etc/profile.d/conda.sh
conda activate snakemake
snakemake --profile prachu_lsf --restart-times 2 --latency-wait 10 #--batch all=1/3
#snakemake -j 4 --restart-times 2 --latency-wait 10
