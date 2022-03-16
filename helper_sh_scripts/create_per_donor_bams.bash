#!/bin/bash


# ${1} is the donor name
# ${2} is the hash_file with donor name and the cell barcodes (with header, tab separated)
# ${3} is temp dir (with sample name as the parent dir)
# ${4} is the output dir (with sample name as the parent dir)
# ${5} is the multiplexed bam file

ml samtools

# For some reason, awk wouln't work
if [ ! -d "${3}" ]; then mkdir -p ${3}; fi
if [ ! -d "${4}" ]; then mkdir -p ${4}; fi

awk -v donor="${1}" '(NR> 1 && $1 == donor){print $2}' ${2} > ${3}${1}.txt
# samtools view -D CB:${3}${1}.txt ${5} -bho "${4}${1}.bam"
samtools view -D CB:${3}${1}.txt ${5} -bho "${3}${1}.bam"
sleep 20
samtools index "${3}${1}.bam" &> /dev/null
sleep 20
samtools view -L "${7}" -o "${4}${1}.bam" "${3}${1}.bam"
sleep 20
samtools index "${4}${1}.bam" &> /dev/null && rm "${3}${1}.bam" "${3}${1}.bam.bai" "${4}${1}.bam.bai" && exit 0 || exit 1
