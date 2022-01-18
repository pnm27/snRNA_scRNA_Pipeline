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
        samtools view -H ${1} > header.sam
	samtools view ${1} | awk -F '\t' -v OFS='\t' -v b="${barc}" '$0~b{print}' > ${1%${base_n}}${l}".sam"
        cat header.sam ${1%${base_n}}${l}".sam" > ${1%${base_n}}${pref}${l}".sam"
        rm ${1%${base_n}}${l}".sam"
        rm header.sam
        samtools view -S -b ${1%${base_n}}${l}".sam" > ${1%${base_n}}${pref}${l}".bam"
        rm ${1%${base_n}}${pref}${l}".sam"

}
done<${2}
date
