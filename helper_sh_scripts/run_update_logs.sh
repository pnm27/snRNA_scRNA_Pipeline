#!/usr/bin/bash
# Runs a python script
# Load necessary Modules


BAM_DIR="/sc/arion/projects/psychAD/STARsolo_bams/"
PICARD_DIR="/sc/arion/projects/psychAD/STARsolo_bams/"
DEMUL_DIR="/sc/arion/projects/psychAD/demultiplex/info/"
MAP_FILE="/sc/arion/projects/psychAD/pnm/Final_out_MAP_2.tsv"
OUTPUT_FILE="/sc/arion/projects/psychAD/pnm/All_logs.tsv"
INPUT_FILE="/sc/arion/projects/psychAD/pnm/fastq_files.txt"
BAM_STRUCT="Sample_<sample>*/"
PICARD_STRUCT="Sample_<sample>*/"
DEMUX_STRUCT="Sample_<sample>*/"


samp_list=""
while read l 
do
{
    # Remove 'cDNA' at the end of the names
    # glob can add cDNA or HTO as required
	samp_name=$(rev <<< ${l} | cut -d "/" -f1 | cut -d "-" -f2- | rev )
	samp_list=${samp_list}" "${samp_name}
}
done < ${INPUT_FILE}

echo "STARTING the script ${0} at:"
date
# Restrict glob expansion
python3 ../helper_py_scripts/update_logs.py -m ${MAP_FILE} -o ${OUTPUT_FILE} -b ${BAM_DIR} -p ${PICARD_DIR} -d ${DEMUL_DIR} ${samp_list} \
 --bam_struct "${BAM_STRUCT}" --pc_struct "${PICARD_STRUCT}" --dem_struct "${DEMUX_STRUCT}" --ss_l --pc_gc --pc_rs --ss_g_f --ss_gf_f --ss_g_s --ss_gf_s --ss_bc --dem_info

echo "FINISHED the script at:"
date