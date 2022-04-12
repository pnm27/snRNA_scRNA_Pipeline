#!/usr/bin/sh
# Runs a python script
# Load necessary Modules


FINAL_COUNT_MATRIX_DIR="/sc/arion/projects/CommonMind/pnm/choroid_plexus/final_count_matrix/solo/"
OUTPUT_DIR="/sc/arion/projects/CommonMind/single-cell-neurogenomics/choroid-plexus/02_h5ad-by-donor/"
SAMPLE_NAME_COL="LAB,SAMPLE,SAMPLE"
LOG_FILE="/sc/arion/projects/CommonMind/pnm/choroid_plexus/All_logs.tsv"
WET_LAB_FILE=""
DONOR_COLUMN="STARsolo,DEMUX,N_CELLS_AFTER_DEMUX_CS"
MULTIPLE_HEADER="3"
DONOR_SEP=","
BATCH_COLUMN="LAB,BATCH,SET"


# IF output_dir doesn't exist, create it
if [ ! -d ${OUTPUT_DIR} ]; then mkdir -p ${OUTPUT_DIR}; fi


if [ -z ${WET_LAB_FILE} ] && [ ! -z ${LOG_FILE} ] && [ ! -z ${DONOR_SEP} ] && [ -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --log_file ${LOG_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP}
	echo "FINISHED the script at:"
	date

elif [ -z ${WET_LAB_FILE} ] && [ ! -z ${LOG_FILE} ] && [ ! -z ${DONOR_SEP} ] && [ ! -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --log_file ${LOG_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP} --batch_col "${BATCH_COLUMN}"
	echo "FINISHED the script at:"
	date

elif [ -z ${WET_LAB_FILE} ] && [ ! -z ${LOG_FILE} ] && [ -z ${DONOR_SEP} ] && [ -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --log_file ${LOG_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER}
	echo "FINISHED the script at:"
	date

elif [ -z ${WET_LAB_FILE} ] && [ ! -z ${LOG_FILE} ] && [ -z ${DONOR_SEP} ] && [ ! -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --log_file ${LOG_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} \
	--batch_col "${BATCH_COLUMN}"
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z ${LOG_FILE} ] && [ ! -z ${DONOR_SEP} ] && [ -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --wet_lab_file ${WET_LAB_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP}
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z ${LOG_FILE} ] && [ ! -z ${DONOR_SEP} ] && [ ! -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --wet_lab_file ${WET_LAB_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} \
	--donor_sep ${DONOR_SEP} --batch_col "${BATCH_COLUMN}"
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z ${LOG_FILE} ] && [ -z ${DONOR_SEP} ] && [ -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --wet_lab_file ${WET_LAB_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} 
	echo "FINISHED the script at:"
	date

elif [ ! -z ${WET_LAB_FILE} ] && [ -z ${LOG_FILE} ] && [ -z ${DONOR_SEP} ] && [ ! -z ${BATCH_COLUMN} ]
then
	echo "STARTING to split the pooled h5ads at:"
	date
	python3 ../helper_py_scripts/split_seth5ad_to_samph5ad.py ${FINAL_COUNT_MATRIX_DIR} ${OUTPUT_DIR} "${SAMPLE_NAME_COL}" --wet_lab_file ${WET_LAB_FILE} --cols "${DONOR_COLUMN}" --multiheader ${MULTIPLE_HEADER} \
	--batch_col "${BATCH_COLUMN}"
	echo "FINISHED the script at:"
	date

else
	echo "None of the conditions required to execute the splitting of the pooled h5ad files is satisfied!"

fi
