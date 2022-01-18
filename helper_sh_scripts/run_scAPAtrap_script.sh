#!/bin/bash

#BSUB -J NPSAD-scAPAtrap
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=5000]
#BSUB -W 48:00
#BSUB -oo /sc/arion/projects/psychAD/pnm/scAPAtrap.stdout
#BSUB -eo /sc/arion/projects/psychAD/pnm/scAPAtrap.stderr
#BSUB -L /bin/bash


ml samtools
ml subread
ml R/4.1.0

sam_path=$(echo $PATH | tr ':' '\n' | grep 'samtools')
fC_path=$(echo $PATH | tr ':' '\n' | grep 'subread')'/featureCounts'
outdir=$(basename ${2})

Rscript /sc/arion/projects/psychAD/pnm/helper_r_scripts/scAPA_script.R ${sam_path} ${fC_path} ${1} ${3} -t ${4} -o ${outdir} -x
