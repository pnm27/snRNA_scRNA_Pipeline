#!/bin/bash


# ${1} is the donor name
# ${2} is the hash_file with donor name and the cell barcodes (with header, tab separated)
# ${3} is temp dir (with sample name as the parent dir) which will contain only barcodes for the respective donors
# ${4} is the output dir (with sample name as the parent dir)
# ${5} is the multiplexed bam file

ml samtools/1.15.1

# To debug:
# 1) Checking inputs
# for val in "$@"; do
#     echo ${val} 
# done
# 2) check all: 
# set -x

# For some reason, awk wouln't work
if [ ! -d "${3}" ]; then mkdir -p ${3}; fi


awk -v donor="${1}" '(NR> 1 && $1 == donor){print $2}' ${2} > ${3}${1}.txt
# samtools view -D CB:${3}${1}.txt ${5} -bho "${4}${1}.bam"
samtools view -D CB:${3}${1}.txt ${5} -bho "${4}${1}.bam" && samtools index "${4}${1}.bam" &> /dev/null && rm "${4}${1}.bam.bai" && exit 0 || exit 1
sleep 20

# Unsetting debug
# set +x
