#!/usr/bin/bash
# ${1} is the bam file (sorted by 'CB' tag)
# ${2} is the barcodes file (called as 'barcodes_list file'), which is tab-separated with donor names as 1st column and corresponding barcodes
# on the 2nd column
# ${2} is assumed to have headers
# ${3} is the dir where the split_bams should go
# ${4} is the temp dir where each file (cell/barcode) is produced with each sample name as a prefix

# example run of this script sh split_bam_indiv_barc_samtools.sh <bam_file> <barcode> <split_bams_dir>
# To do: loop over the donors (1st column; use awk's arrays) and get the corresponding barcodes
# Default naming style of split bams: <split_bams_dir><sample_name>/<donor>.bam 
# sample name will be the basename extracted from the bam file
ml samtools

sample_name=$(basename ${1} | sed 's/_Aligned.sortedByCB.out.bam//g')
bam_dir=$(dirname ${1})
echo "Splitting bam by barcodes at: "
date

# Create split_bams dir per donor to facilitate parallelization
split_bam_dir=$(op=("${3}" "${sample_name}"); printf '%s/' "${op[@]%/}")
if [ ! -d "${split_bam_dir}" ]; then echo "Entered directory to store split bams doesn't exists! Creating it."; mkdir -p ${split_bam_dir}; fi
# Make temp dir to facilitate parallelization
if [ ! -d "${4}" ]; then echo "Entered temp directory doesn't exists! Creating it."; mkdir -p ${4}; fi
# Make temp dir per donor to facilitate parallelization

# Temp file with 3 columns containing barcode, line number of first read, line number of last read

# Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
split_n=800

# Skip these lines for each search
# For the first run skip headers
headers_sam=$(echo "${4%/}"/"${sample_name}"_"headers.sam")
samtools view -H ${1} > ${headers_sam}

# To split bam by barcodes
samtools view ${1} | while read l
do
{
	# Stores current value of 'CB' tag
	bc=$(echo "${l}" | sed -nr "s/.*CB:Z:([ATGCN]+).*/\1/p")

	# If no "CB" tag is found then skip
	if [ -z "${bc}" ]
	then
		continue
	fi

	temp_b="${bc}.bam"
	temp_s="${bc}.sam"
	t_op_bam=$(echo "${4%/}"/"${sample_name}"_"${temp_b}")
	t_op_sam=$(echo "${4%/}"/"${sample_name}"_"${temp_s}")
	# If barcode is already found (variable called 'old_bc' will be equal to the variable 'bc')
	# and read is not already present in the sam file

	# if 'bc' is not in barcodes_list file and 'old_bc is unset' then skip the read
	if [ -z "${old_bc}" ] && ! grep -q "${bc}" "${2}"
	then
		continue

	# If both 'old_bc' (can't be empty or unset) and 'bc' are the same and aren't in the barcodes_list file then skip the read
	elif [ ! -z "${old_bc}" ] && [ "${bc}" == "${old_bc}" ] && ! grep -q "${bc}" "${2}"
	then
		continue

	# If barcode is already found i.e. first read for this barcode was already found by this loop, previously
	# And the sam exists for that barcode (if the barcode is already found then it should have sam-the build of this loop always guarantees
	# the formation of sam file for the "filtered" barcode i.e the bc is present in the barcodes_list file)
	elif [ ! -z "${old_bc}" ] && [ "${old_bc}" == "${bc}" ] && [ -f "${t_op_sam}" ]
	then
	{
		# echo "Old bc is: ${old_bc}"
		# echo "current bc is: ${bc}"
		# echo "bc_found loop"
		# echo "${l}"
		echo "${l}" >> ${t_op_sam}
	}

	# when "new" barcode is found (old_bc can't be empty/unset), 'old_bc' is in the barcodes_list file, and 'bc' doesn't exist
	# in the barcodes_list file then create the bam for the 'old_bc', set 'old_bc' to 'bc'
	elif [ ! -z "${old_bc}" ] && [ "${bc}" != "${old_bc}"  ] && grep -q "${old_bc}" "${2}" && ! grep -q "${bc}" "${2}"
	then
	{
		# echo "Old bc is: ${old_bc}"
		# echo "current bc is: ${bc}"
		prev_bam=$(echo "${4%/}"/"${sample_name}"_"${old_bc}.bam")
		prev_sam=$(echo "${4%/}"/"${sample_name}"_"${old_bc}.sam")
		# echo "new_bc loop"
		# head "${prev_sam}"
		# tail "${prev_sam}"
		# cat "${headers_sam}" "${prev_sam}" > "${old_bc}_check.sam"
		samtools view -@ 6 -bS <(cat "${headers_sam}" "${prev_sam}") > ${prev_bam} # 'S' option is for being "safe"-confirm the input is sam
		rm ${prev_sam}
		
	}


	# If bc not found before by this script (either first barcode present in the barcodes_list file or
	# "new" barcode wrt loop inside bam file)	
	elif grep -q "${bc}" "${2}"
	then
	{
		# when "new" barcode is found (old_bc can't be empty/unset), 'old_bc' is in the barcodes_list file, and 'bc' exists
		# in the barcodes_list file then create the bam for the 'old_bc' 
		if [ ! -z "${old_bc}" ] && [ "${old_bc}" != "${bc}" ] && grep -q "${old_bc}" "${2}"
		then
		{
			# echo "Old bc is: ${old_bc}"
			# echo "current bc is: ${bc}"
			prev_bam=$(echo "${4%/}"/"${sample_name}"_"${old_bc}.bam")
			prev_sam=$(echo "${4%/}"/"${sample_name}"_"${old_bc}.sam")
			# echo "new_bc loop"
			# head "${prev_sam}"
			# tail "${prev_sam}"
			# cat "${headers_sam}" "${prev_sam}" > "${old_bc}_check.sam"
			samtools view -@ 6 -bS <(cat "${headers_sam}" "${prev_sam}") > ${prev_bam} # 'S' option is for being "safe"-confirm the input is sam
			rm ${prev_sam}
		}
		fi

		# If sam is present and read not present
		if [ -f "${t_op_sam}" ] && ! grep -Fxoq "${l}" "${t_op_sam}"
		then
			echo "${l}" >> ${t_op_sam}

		# If sam is present and read present
		elif [ -f "${t_op_sam}" ] && grep -Fxoq "${l}" "${t_op_sam}"
		then
			continue

		# if sam not present but bam is present and read present
		elif [ ! -f "${t_op_sam}" ] && [ -f "${t_op_bam}" ] && samtools view "${t_op_bam}" | grep -Fxoq "${l}"
		then
			continue

		# If sam is not present but bam is present and read not present	
		elif [ ! -f "${t_op_sam}" ] && [ -f "${t_op_bam}" ] && ! samtools view "${t_op_bam}" | grep -Fxoq "${l}"
		then
			echo "${l}" | cat <(samtools view ${t_op_bam}) -  > ${t_op_sam}

		#  If sam and bam both not present
		elif [ ! -f "${t_op_sam}" ] && [ ! -f "${t_op_bam}" ]
		then
			echo "${l}" > ${t_op_sam}

		# Not expecting any new cases
		else
			continue
			# echo "Wrong action!"
		fi

	}
	else
	{
		continue
		# echo "Old bc is: ${old_bc}"
		# echo "current bc is: ${bc}"
		# echo "missed some condition"			
	}
	fi

	old_bc=${bc}

}
done
echo "Finished splitting by barcodes at: "
date

# Remove all intermediate sam files
sam_files=$(echo "${4%/}"/"*.sam")
rm ${sam_files}

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
			samtools merge -@ 8 -o ${temp_bam} ${inp_files}
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

# Remove bams per barcode for the sample (only retain donor-level)
# to_rem="${4%/}/${sample_name}_*.bam"
# eval rm "${to_rem}"

echo -n "Finished Processing "${sample_name}" at: "
date
