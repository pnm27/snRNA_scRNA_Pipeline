#!/usr/bin/bash
# ${1} is the bam file
# ${2} is assumed to have headers
# ${2} is the barcodes file, which is tab-separated with donor names as 1st column and corresponding barcodes
# on the 2nd column
# ${3} is the dir where the split_bams should go
# ${4} is the temp dir

# example run of this script sh split_bam_indiv_barc_samtools.sh <bam_file> <barcode> <split_bams_dir>
# To do: loop over the donors (1st column; use awk's arrays) and get the corresponding barcodes
# Default naming style of split bams: <split_bams_dir><sample_name>/<donor>.bam 
# sample name will be the basename extracted from the bam file
ml samtools

sample_name=$(basename ${1} | sed 's/_Aligned.sortedByCoord.out.bam//g')
date

# Create split_bams dir per donor to facilitate parallelization
temp_dir=$(op=("${3}" "${sample_name}"); printf '%s/' "${op[@]%/}")
if [ ! -d "${temp_dir}" ]; then echo "Entered temp directory doesn't exists! Creating it."; mkdir -p ${temp_dir}; fi
# Make temp dir to facilitate parallelization
if [ ! -d "${4}" ]; then echo "Entered temp directory doesn't exists! Creating it."; mkdir -p ${4}; fi
# Make temp dir per donor to facilitate parallelization



# Read as array
awk 'NR>1{a[$1]++}END{for (b in a) print b}' ${2} | while read l 
do
{
	echo ${l}
	awk -v donor="${l}" '(NR>1 && $1 == donor){print $2}' ${2} | while read bc
	do
	{
		# pref=${3}
		# barc="CB:Z:${bc}"
		temp_f=temp_${bc}.bam
		# t_op_fi=$(op=("${4}" "${temp_f}"); printf '%s/' "${op[@]%/}" | sed "s/\/$//" )
		t_op_fi=$( echo "${4%/}"/"${sample_name}"_"${temp_f}" )
		samtools view -bhd "CB:${bc}" ${1} > ${t_op_fi}
	}
	done
	t_bams="${4%/}/${sample_name}*.bam"
	inp_files=$( ls ${t_bams} | tr '\n' ' ' )
	op_f=$( echo "${temp_dir%/}" | sed "s/$/\/${l}\.bam/"  )
	samtools merge -o ${op_f} ${inp_files}

	# Clear temp files
	rm ${inp_files}
}
done
echo -n "Finished Processing "${sample_name}" at: "
date
