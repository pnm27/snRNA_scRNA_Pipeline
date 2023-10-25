#! /usr/bin/sh

#BSUB -J AMP_PD_agg_mat
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=8000]
#BSUB -W 08:00
#BSUB -oo /sc/arion/projects/CommonMind/pnm/AMP_PD/test_new_script.stdout
#BSUB -eo /sc/arion/projects/CommonMind/pnm/AMP_PD/test_new_script.stderr
#BSUB -L /bin/bash


ml anaconda3
ml R/4.1.0

source /sc/arion/work/prashf01/conda/envs/pegasus/etc/profile.d/conda.sh
conda activate pegasus

python3 create_aggregate_matrix.py /sc/arion/projects/CommonMind/pnm/AMP_PD/final_count_matrix/test/ \
 /sc/arion/projects/psychAD/pnm/Hs_allchr_MT.txt /sc/arion/projects/CommonMind/pnm/AMP_PD/test_combo.h5ad \
 --conversion_file conversion_vir_w_gt2.tsv --conversion_file_cols Sample,Donor_name,donor,Final_breg \
 --suffix w_gt --h5ad_col donor_,brain_reg_ --conversion_file conversion_vir_first_run2.tsv \
 --conversion_file_cols Sample,Donor_name,donor,Final_breg --suffix wo_gt --h5ad_col donor_,brain_reg_ \
 --column_prefix SubID_vs_ -g GRCh38
