#!/usr/bin/bash
# ${1} is the bam file
# ${2} is assumed to have headers
# ${2} is the barcodes file, which is tab-separated with donor names as 1st column and corresponding barcodes
# on the 2nd column

# To do: loop over the donors (1st column; use awk's arrays) and get the corresponding barcodes
ml samtools

sample_name=$(basename ${1} | sed 's/_Aligned.sortedByCoord.out.bam//g')
date
awk -F "\t" 'NR>1{a[$1]++}' ${2} | while read l # Read as array
do
{
	awk -F "\t" -v donor="${l}" '(NR>1 && $1 == donor){print $2}' ${2} | while read bc
	do
	{
		pref=${3}
		barc="CB:Z:${bc}"
		samtools view -bhd "CB:"${bc} ${1} > ${pref}${l}".bam"
	}
	done
}
done
echo -n "Finished Processing "${sample_name}" at: "
date
