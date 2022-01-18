#!/bin/bash

#BSUB -J NPSAD-new_split_bam
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=5000]
#BSUB -W 48:00
#BSUB -oo /sc/arion/projects/psychAD/pnm/new_split_bam.stdout
#BSUB -eo /sc/arion/projects/psychAD/pnm/new_split_bam.stderr
#BSUB -L /bin/bash
