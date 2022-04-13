#!/usr/bin/bash

FASTQS_DIR="/sc/arion/projects/CommonMind/pnm/choroid_plexus/fastq_links/"
CHECKSUM_DIR="checksum/"
FASTQS_SUFFIX="fastq.gz"
LOG_FILE="/sc/arion/projects/CommonMind/pnm/choroid_plexus/hash_check_fastqs.txt"

# echo "This script is associated with the snakemake pipeline. Hence, it assumes the files are links"
# echo "This script assumes that cDNA and HTO files differ only in the term (cDNA or HTO) in their path."
# echo "e.g. if the path to R1 fastq of the cDNA is fastq_dir/sample_1_cDNA/sample_1_cDNA.R1.fastq.gz then"
# echo "path to R1 fastq of the HTO is assumed to be fastq_dir/sample_1_HTO/sample_1_HTO.R1.fastq.gz"
# echo
echo "By default, this script assumes the checksum_dir to be present in the same dir as the fastqs"
echo "and the hashes of the file have same basename but appended with the suffex '.md5'"
echo "e.g. if the path to R1 fastq of the cDNA is fastq_dir/sample_1_cDNA/sample_1_cDNA.R1.fastq.gz and checksum_dir is 'checksum/' then"
echo "path to the md5 file of the R1 fastq of the cDNA is assumed to be fastq_dir/sample_1_HTO/checksum/sample_1_HTO.R1.fastq.gz.md5"
echo
echo "Only if a file's hash doesn't match it will be logged into the log file"

# Touch the log_file, if it doesn't exists (similarly for its parent folder)
if [ ! -d "${LOG_FILE%$(basename ${LOG_FILE})}" ]; then mkdir -p "${LOG_FILE%$(basename ${LOG_FILE})}"; fi
if [ ! -s "${LOG_FILE}" ]; then touch "${LOG_FILE}"; fi


find -L "${FASTQS_DIR}" -name "*${FASTQS_SUFFIX}" -exec readlink -f {} \; | while read l
do
	b=$(basename ${l})
	checksum_file="${l%${b}}${CHECKSUM_DIR%/}/${b}.md5"
	if [ ! -f "${checksum_file}" ]; then
		echo "Checksum file couldn't be found!" >> ${LOG_FILE}
		echo "File: ${l}" >> ${LOG_FILE}
		echo "Checksum File: ${checksum_file}" >> ${LOG_FILE}
		echo "Hence, skipping this file" >> ${LOG_FILE}
		continue
	fi

	hash_file_value=$(cat ${checksum_file} | cut -d " " -f1)
	hash_value=$(md5sum ${l} | cut -d " " -f1)
	if [[ "${hash_value}" != "${hash_file_value}" ]]; then
		echo "Hash_values for the file ${b} weren't the same with the one (${checksum_file}) in the checksum dir (${CHECKSUM_DIR})" >> ${LOG_FILE}
	fi

done

# If the log file is still of size '0' then no problems occurred!
if [ ! -s "${LOG_FILE}" ]; then echo "All went well!" >> "${LOG_FILE}"; fi