#!/usr/bin/bash

###
#If the script is executed as:
# sh create_inp_for_snkmk.sh 1 5 /sc/arion/projects/psychAD/pnm/fastq_files.txt
#Then remove existing fastq_files.txt and  create fastq_files.txt containing all samples from round1 to round5 for snakemake's run
###

start_round=${1}
end_round=${2}

if [ -f "${3}" ]; then rm ${3}; fi

seq ${start_round} ${end_round} | while read l
do
{
	ls /sc/arion/projects/psychAD/fastqs/round${l}/*-cDNA/*.R1.fastq.gz | sed "s/^\/sc\/arion\/projects\/psychAD\/fastqs\///g" | sed "s/.R1.fastq.gz//g" >> /sc/arion/projects/psychAD/pnm/fastq_files.txt

}
done
