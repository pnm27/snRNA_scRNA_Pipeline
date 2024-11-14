#!/usr/bin/bash
# Runs a python script
# Load necessary Modules

declare -A params # contains the parameter for python script
declare -A args # contains the values for above params

# DON'T CHANGE THE BELOW SET
params["BAM_DIR"]="-b"
params["PICARD_DIR"]="-p"
params["DEMUL_DIR"]="-d"
params["MAP_FILE"]="-m"
params["OUTPUT_FILE"]="-o"
params["BAM_STRUCT"]="--bam_struct"
params["PICARD_STRUCT"]="--pc_struct"
params["DEMUX_STRUCT"]="--dem_struct"

# CHANGE THIS ACCORDING TO PROJECT REQUIREMENTS
INPUT_FILE="/sc/arion/projects/psychAD/pnm/fastq_files.txt"
args["BAM_DIR"]="/sc/arion/projects/psychAD/STARsolo_bams/" # REQUIRED
args["PICARD_DIR"]="/sc/arion/projects/psychAD/STARsolo_bams/"
args["DEMUL_DIR"]="/sc/arion/projects/psychAD/demultiplex/info/"
args["MAP_FILE"]="/sc/arion/projects/psychAD/pnm/Final_out_MAP_2.tsv"
args["OUTPUT_FILE"]="/sc/arion/projects/psychAD/pnm/All_logs.tsv" # REQUIRED
args["BAM_STRUCT"]="Sample_<sample>*/" # REQUIRED
args["PICARD_STRUCT"]="Sample_<sample>*/"
args["DEMUX_STRUCT"]="Sample_<sample>*/"


samp_list=""
while read l 
do
{
    if ! grep -q "^#" <<< ${l}; then
        # Remove 'cDNA' at the end of the names
        # glob can add cDNA or HTO as required
        samp_name=$(rev <<< ${l} | cut -d "/" -f1 | cut -d "-" -f2- | rev )
        samp_list=${samp_list}" "${samp_name}
    fi
}
done < ${INPUT_FILE}

# build command
cmd_str=""
for key in "${!args[@]}"
do
    cmd_str+=$([ ! -z "${args[$key]}" ] && echo " ${params[$key]} ${args[$key]}")
done

echo "STARTING the script ${0} at:"
date
echo
set -x
# Restrict glob expansion
python3 ../helper_py_scripts/update_logs.py ${cmd_str} ${samp_list} \
    --ss_l --ss_g_f --ss_gf_f --ss_g_s --ss_gf_s --ss_bc --pc_gc --pc_rs --dem_info
    
set +x
echo
echo "FINISHED the script at:"
date