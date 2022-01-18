#!/usr/bin/bash
# ${1} is bam file
# ${2} is the barcodes file (in the future we can have barcode associated with each donor)


ml samtools

date
while read l
do
{
	base_n=$( basename ${1})
	pref=$( echo ${base_n} | sed 's/_L.*.bam/_/g' )   # Just to be safe 'L' doesn't occur anywhere else
	barc="CB:Z:${l}"
	samtools view -bhd "CB:"${l} ${1} > ${1%${base_n}}${pref}${l}"_samtools.bam"

}
done<${2}
date
