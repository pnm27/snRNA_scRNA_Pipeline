#!/usr/bin/bash
# Input:
# ${1} is the bam file (sorted by 'CB' tag)

# Output:
# ${2} is the file containing barcode and number of reads per barcode 
# ${3} is the shell script which will split the bam
# ${4} is the split bams dir

echo "Starting script at: "
date

ml samtools

# For setting folder structure
sample_name=$(basename ${1} | sed 's/_Aligned.sortedByCB.out.bam//g')
bam_dir=$(dirname ${1})
round_dir=$(basename $(dirname ${bam_dir}))
sam_f="${1::-3}sam"


echo "Sample name is: ${sample_name}"
echo "Bam dir is: ${bam_dir}"
echo "round num is: ${round_dir}"

# Create split_bams dir per sample to facilitate parallelization
# split_bam_dir=$(op=("${4}" "${round_dir}" "${sample_name}"); printf '%s/' "${op[@]%/}")
split_bam_dir=${4}
if [ ! -d "${split_bam_dir}" ]; then echo "Entered directory to store split bams doesn't exists! Creating it."; mkdir -p ${split_bam_dir}; fi

# Sam and the barcodes files don't exist
if [ ! -f "${sam_f}" ] && [ ! -f "${2}" ]
then
{
    echo "Sam and barcodes files are absent! Creating sam and creating barcodes file at: "
    date

    samtools view ${1} > "${sam_f}"

    echo "Sam Created at: "
    date
    grep -oE "CB:Z:([ATGCN])+" "${sam_f}" | sort | uniq -c > ${2} && echo "Barcodes file at: " && date
}
# Sam exists but the barcodes file don't exist
elif [ -f "${sam_f}" ] && [ ! -f "${2}" ]
then
{
    echo "Sam file present but barcodes file is not! Creating Barcodes file at: "
    date
    grep -oE "CB:Z:([ATGCN])+" "${sam_f}" | sort | uniq -c > ${2} && echo "Barcodes file at: " && date
}
# Sam exists and the barcodes file exist
elif [ -f "${sam_f}" ] && [ -f "${2}" ]
then
{
    echo "Sam file and Barcodes file both exists."
    date
}
# Sam file doesn't exist but barcodes file exists
else
{
    echo "Sam file doesn't exist but barcodes file does! Creating sam but not the barcodes file at: "
    date

    samtools view ${1} > "${sam_f}"

    echo "Sam Created at: "
    date   
}
fi


(
cat<<EOF
#!/usr/bin/bash

ml samtools

echo "Starting to split bams at: "
date

if [ ! -f "${bam_dir}/${sample_name}_headers.sam" ]; then samtools view -H ${1} > ${bam_dir}/${sample_name}_headers.sam; fi

(
EOF
) > ${3}


# echo -e "samtools view -H ${1} > ${bam_dir}/${sample_name}_headers.sam\n\n(" >> ${3}


while read -a l
do
{
    n_reads=${l[0]}
    bc=${l[1]:5} # The barcodes are captured as "CB:Z:[ATGCN]+"
    echo "head -n ${n_reads} | cat "${bam_dir}/${sample_name}_headers.sam" - | samtools view -bS - > ${split_bam_dir}${bc}.bam" >> ${3}

}
done<${2}

echo -e ")<${sam_f}\n\nrm ${bam_dir}/${sample_name}_headers.sam\n\nrm ${sam_f}\n\necho "Finished splitting bams at: "\ndate" >> ${3}
# echo -e ")<${sam_f}\n\necho "Finished splitting bams at: "\ndate" >> ${3}

echo "Script finished at: "
date
