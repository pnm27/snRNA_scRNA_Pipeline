#!/usr/bin/bash
# Runs a python script
# Load necessary Modules


FINAL_COUNT_MATRIX_DIR="/sc/arion/projects/CommonMind/pnm/choroid_plexus/final_count_matrix/solo/"
OUTPUT_DIR="/sc/arion/projects/CommonMind/single-cell-neurogenomics/choroid-plexus/02_h5ad-by-donor"
SAMPLE_NAME_COL="LAB SAMPLE SAMPLE"
LOG_FILE="/sc/arion/projects/CommonMind/pnm/choroid_plexus/All_logs.tsv"
WET_LAB_FILE=""
DONOR_COLUMN="STARsolo DEMUX N_CELLS_AFTER_DEMUX_CS"
MULTIPLE_HEADER="3"
DONOR_SEP=","
SAMPLE_NAME_FMT="-([0-9]+-[A-Za-z0-9]+)+"
SAMPLE_NAME_GROUPS="1"


if [ -z ${WET_LAB_FILE} ] && [ ! -z $ {LOG_FILE} ] && [ ! -z ${SAMPLE_NAME_FMT} ] && [ ! -z ${SAMPLE_NAME_GROUPS} ] && [ ! -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --log_file ${LOG_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP} --sample_name_format ${SAMPLE_NAME_FMT} --samp_name_grp ${SAMPLE_NAME_GROUPS}
	echo "FINISHED the script at:"
	date

elif [ -z ${WET_LAB_FILE} ] && [ ! -z $ {LOG_FILE} ] && [ -z ${SAMPLE_NAME_FMT} ] && [ ! -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --log_file ${LOG_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP}
	echo "FINISHED the script at:"
	date

elif [ -z ${WET_LAB_FILE} ] && [ ! -z $ {LOG_FILE} ] && [ ! -z ${SAMPLE_NAME_FMT} ] && [ ! -z ${SAMPLE_NAME_GROUPS} ] && [ -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --log_file ${LOG_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER} \
	--sample_name_format ${SAMPLE_NAME_FMT} --samp_name_grp ${SAMPLE_NAME_GROUPS}
	echo "FINISHED the script at:"
	date

elif [ -z ${WET_LAB_FILE} ] && [ ! -z $ {LOG_FILE} ] && [ -z ${SAMPLE_NAME_FMT} ] && [ -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --log_file ${LOG_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER}
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z $ {LOG_FILE} ] && [ ! -z ${SAMPLE_NAME_FMT} ] && [ ! -z ${SAMPLE_NAME_GROUPS} ] && [ ! -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --wet_lab_file ${WET_LAB_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP} --sample_name_format ${SAMPLE_NAME_FMT} --samp_name_grp ${SAMPLE_NAME_GROUPS}
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z $ {LOG_FILE} ] && [ -z ${SAMPLE_NAME_FMT} ] && [ ! -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --wet_lab_file ${WET_LAB_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP}
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z $ {LOG_FILE} ] && [ ! -z ${SAMPLE_NAME_FMT} ] && [ ! -z ${SAMPLE_NAME_GROUPS} ] && [ -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --wet_lab_file ${WET_LAB_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER} \
	--sample_name_format ${SAMPLE_NAME_FMT} --samp_name_grp ${SAMPLE_NAME_GROUPS}
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z $ {LOG_FILE} ] && [ -z ${SAMPLE_NAME_FMT} ] && [ -z ${DONOR_SEP} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} ${SAMPLE_NAME_COL} --wet_lab_file ${WET_LAB_FILE} --cols ${DONOR_COLUMN} --multiheader ${MULTIPLE_HEADER}
	echo "FINISHED the script at:"
	date
else
	echo "None of the conditions required to execute the splitting of the pooled h5ad files is satisfied!"

fi
