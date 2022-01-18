#!/usr/bin/bash
# ${1} is the barcodes file (called as 'barcodes_list file'), which is tab-separated with donor names as 1st column and corresponding barcodes
# on the 2nd column
# ${1} is assumed to have headers
# ${2} is the dir where the split_bams should go
# ${3} is the dir where each per-cell bams are there
# ${4} is the numerical limit to split the barcodes file (to avoid opening a lot of bam files)

# This script loops over the donors (1st column; using awk's arrays) and get the corresponding barcodes
# Default naming style of split bams: <split_bams_dir><donor>.bam . NOTE: <split_bams_dir> has a file structure associated with it
# e.g.(/sc/arion/projects/psychAD/STARsolo_split_filt_bams/round1/NPSAD-20201013-A1-cDNA/M34966.bam)
#      ----------------------split_bams_dir--------------------------------------------- -donor-
# sample name will be the basename extracted from the bam file

# Split barcodes if more than this number belonging to the same donor (can't merge files more than what specified by `ulimit -n`)
echo "Starting to assemble bams at donor-level at: "
date

split_n=${4}

split_bam_dir=${2}

if [ ! -d "${split_bam_dir}" ]; then mkdir -p "${split_bam_dir}"; fi

ml samtools

# To assign bc to each donor
awk 'NR>1{a[$1]++}END{for (b in a) print b}' ${1} | while read d
do
{
	n_bcs=$(awk -v donor="${d}" '(NR>1 && $1 == donor){a++}END{print a}' ${1})
	per_bc_bams_dir=$(echo "${3%/}/") # Directory where bams per cell(barcode) is present
	op_f=$( echo "${split_bam_dir%/}" | sed "s/$/\/${d}\.bam/"  )
	if [ "${n_bcs}" -lt "$((split_n+1))" ]
	then
	{
		# inp_files=$( ls ${t_bams} | tr '\n' ' ' )
		inp_files=$(awk -v donor="${d}" -v pd="${per_bc_bams_dir}" 'BEGIN{ORS=" ";}(NR>1 && $1 == donor){print pd$2".bam"}' ${1})
		
		samtools merge -@ 8 -o ${op_f} ${inp_files}

		# Clear all per-bc bams
		rm ${inp_files}
	}
	else
	{
		# Temp txt file containing all bc per donor
		t_op_t=$(echo "${split_bam_dir%/}" | sed "s/$/\/${d}_bc\.txt/")
		# Prefix for each file split file (split the file that contains all bc per donor)
		t_op_spl=$(echo "${split_bam_dir%/}" | sed "s/$/\/${d}_bc_/")
		awk -v donor="${d}" '(NR>1 && $1 == donor){print $2}' ${1} > ${t_op_t}
		split -l ${split_n} ${t_op_t} ${t_op_spl}
		# List of all split files
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
				inp_files=$(awk -v pd="${per_bc_bams_dir}" -v m_bam="${temp_bam}" 'BEGIN{ORS=" ";}{print pd$1".bam"}END{print m_bam}' ${file})
				c+=1
			}
			else
				inp_files=$(awk -v pd="${per_bc_bams_dir}" 'BEGIN{ORS=" ";}{print pd$1".bam"}' ${file})
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