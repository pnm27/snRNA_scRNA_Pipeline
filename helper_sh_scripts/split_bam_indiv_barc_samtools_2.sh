#!/usr/bin/bash
# ${1} is the bam file
# ${2} is assumed to have headers
# ${2} is the barcodes file, which is tab-separated with donor names as 1st column and corresponding barcodes
# on the 2nd column
# ${3} is the dir where the split_bams should go
# ${4} is the temp dir where each file (cell/barcode) is produced inside each sample (folder)

# example run of this script sh split_bam_indiv_barc_samtools.sh <bam_file> <barcode> <split_bams_dir>
# To do: loop over the donors (1st column; use awk's arrays) and get the corresponding barcodes
# Default naming style of split bams: <split_bams_dir><sample_name>/<donor>.bam 
# sample name will be the basename extracted from the bam file
ml samtools

sample_name=$(basename ${1} | sed 's/_Aligned.sortedByCB.out.bam//g')
bam_dir=$(dirname ${1})
date

# Create split_bams dir per donor to facilitate parallelization
split_bam_dir=$(op=("${3}" "${sample_name}"); printf '%s/' "${op[@]%/}")
if [ ! -d "${split_bam_dir}" ]; then echo "Entered directory to store split bams doesn't exists! Creating it."; mkdir -p ${split_bam_dir}; fi
# Make temp dir to facilitate parallelization
if [ ! -d "${4}" ]; then echo "Entered temp directory doesn't exists! Creating it."; mkdir -p ${4}; fi
# Make temp dir per donor to facilitate parallelization

# Temp file with 3 columns containing barcode, line number of first read, line number of last read

# Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
split_n=2


# Skip these lines for each search
# For the first run skip headers
n_header=$(samtools view -H ${1} | wc -l)
declare -i skip_rows=0

# To split bam by barcodes
tail -n +2 ${2} | tr '\t' ' ' | cut -f2 | while read -a bc
do
{
	temp_f="${bc[1]}.bam"
	t_op_bam=$(echo "${4%/}"/"${sample_name}"_"${temp_f}")
	# echo ${skip_rows}
	# echo ${bc[1]}
	# echo ${t_op_bam}
	if [ "${skip_rows}" -eq 0 ]
	then
		first_occur=$(samtools view ${1} | grep -n -m1 "CB:Z:${bc[1]}" | cut -d ":" -f1)
	else
		first_occur=$(samtools view ${1} | sed "1,${skip_rows}d" - | grep -n -m1 "CB:Z:${bc[1]}" | cut -d ":" -f1)
	fi

	samtools view -h ${1} | grep -E "^@|CB:Z:${bc[1]}" | samtools view -b - > ${t_op_bam} && bc_reads=$(samtools view ${t_op_bam} | wc -l)
	skip_rows=$((skip_rows+first_occur+bc_reads-1))
	# echo ${skip_rows}
	# echo -e "\n\n"
}
done
echo "Finished splitting by barcodes"
date



# Donors with more than 800 bc assigned to them:
# sed -e '/Subj_ID/d' <(cut -f1 ${2} | sort | uniq -c) | awk '$1>800{print $2}'


# To assign bc to each donor
awk 'NR>1{a[$1]++}END{for (b in a) print b}' ${2} | while read d
do
{
	n_bcs=$(awk -v donor="${d}" '(NR>1 && $1 == donor){a++}END{print a}' ${2})
	pref_bams="${4%/}/${sample_name}_"
	op_f=$( echo "${split_bam_dir%/}" | sed "s/$/\/${d}\.bam/"  )
	if [ "${n_bcs}" -lt "$((split_n+1))" ]
	then
	{
		# inp_files=$( ls ${t_bams} | tr '\n' ' ' )
		inp_files=$(awk -v donor="${d}" -v pd="${pref_bams}" 'BEGIN{ORS=" ";}(NR>1 && $1 == donor){print pd$2".bam"}' ${2})
		
		samtools merge -@ 8 -o ${op_f} ${inp_files}

		# Clear temp files
		rm ${inp_files}
	}
	else
	{
		t_op_t=$(echo "${split_bam_dir%/}" | sed "s/$/\/${d}_bc\.txt/")
		t_op_spl=$(echo "${split_bam_dir%/}" | sed "s/$/\/${d}_bc_/")
		awk -v donor="${d}" '(NR>1 && $1 == donor){print $2}' ${2} > ${t_op_t}
		split -l ${split_n} ${t_op_t} ${t_op_spl}
		spl_f="${t_op_spl}*"

		# Counter for intermediate "merged_files"
		declare -i c=1

		for file in ${spl_f}
		do
		{

			if [ "${c}" -gt 1 ]
			then
			{
				c=c-1
				temp_bam=$(echo "${split_bam_dir%/}" | sed "s/$/\/temp_${d}_${c}\.bam/")
				inp_files=$(awk -v pd="${pref_bams}" -v m_bam="${temp_bam}" 'BEGIN{ORS=" ";}{print pd$1".bam"}END{print m_bam}' ${file})
				c+=1
				# echo ${temp_bam}
				# echo ${inp_files}
			}
			else
				inp_files=$(awk -v pd="${pref_bams}" 'BEGIN{ORS=" ";}{print pd$1".bam"}' ${file})
			fi
		
			temp_bam=$(echo "${split_bam_dir%/}" | sed "s/$/\/temp_${d}_${c}\.bam/")
			samtools -@ 8 merge -o ${temp_bam} ${inp_files}
			c+=1

		}
		done
		# Remove temp files
		rm ${spl_f}
		rm ${t_op_t}


		# Rename last "merged" bam file as the final bam file
		mv ${temp_bam} ${op_f}

		# Remove all temp "merged" bams
		to_rem=$(echo "${split_bam_dir%/}" | sed "s/$/\/temp_${d}_/")
		eval rm "${to_rem}*.bam"

	}
	fi

}
done

# Remove all split bams for the sample
to_rem="${4%/}/${sample_name}_*.bam"
eval rm "${to_rem}"

echo -n "Finished Processing "${sample_name}" at: "
date
