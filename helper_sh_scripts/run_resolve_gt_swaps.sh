#!/usr/bin/sh


#BSUB -J run_resolve_gt_swaps_NPSAD
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 2
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=4000]
#BSUB -W 48:00
#BSUB -oo /sc/arion/projects/psychAD/pnm/resolve_gt_swaps.stdout
#BSUB -eo /sc/arion/projects/psychAD/pnm/resolve_gt_swaps.stderr
#BSUB -L /bin/bash


ml purge
ml anaconda3
ml R/4.1.0

python3 /sc/arion/work/prashf01/NPSAD_test3/helper_py_scripts/resolve_gt_swaps.py