#!/usr/bin/bash

###
#Create input files called r1_samples, r2_samples, etc
#eg to create r5_samples:
# cd /sc/arion/projects/psychAD/Single_cell_data/fastq/round5
# find -mindepth 1 -maxdepth 1 -type d -printf "%P\n" > ../../../r5_samples

#Use the files called r1_samples, r2_samples..., which contains names like "Sample_NPSAD-20201105-C2-cDNA", "Sample_NPSAD-20201022-A2-HTO", etc., as one of the inputs
#For each rounds, go into the respective round's folder (the one inside /sc/arion/projects/psychAD/fastqs ) and run this

#Example run (for samples in "round1" folder):
# sh create_r1_samp_links.sh r1_sample round1
###


sed 's/Sample_//' ${1} | while read l
do
{
   # This variable contains something like this /sc/arion/projects/psychAD/Single_cell_data/fastq/round1/Sample_NPSAD-20201103-C2-cDNA/fastq/NPSAD-20201103-C2-cDNA_AACCACGCAT_HNYM3DSXY_L001_001
   fpath=$(echo $(du -h "/sc/arion/projects/psychAD/Single_cell_data/fastq/"${2}"/Sample_"${l}/fastq/${l}*R1.fastq.gz | sort -rh | head -1 | cut -f2) | sed "s/.R1.fastq.gz//g")
   # This variable will then contain NPSAD-20201103-C2-cDNA_AACCACGCAT_HNYM3DSXY_L001_001 or NPSAD-20201125-A2-cDNA_CTCTAGCGAG-GATGAAGATA_HVVKTDSXY_L001_001
   temp=$( basename ${fpath})
   # Assuming the naming is consistent, this will name the link as NPSAD-20201103-C2-cDNA_L001_001
   fn=$( cut -d "_" -f1 <<< ${temp})"_"$( cut -d "_" -f4 <<< ${temp})"_"$( cut -d "_" -f5 <<< ${temp})

   # For each sample create its own folder and if the rounds folder is not present then create it too
   if [ ! -d "/sc/arion/projects/psychAD/fastqs/${2}/Sample_${l}" ]; then mkdir -p "/sc/arion/projects/psychAD/fastqs/${2}/Sample_${l}" ; fi

   # Assumption is that the same basename is correct for both R1 and R2 files (logically it should be)
   ln -s ${fpath}".R1.fastq.gz" "/sc/arion/projects/psychAD/fastqs/"${2}"/Sample_"${l}/${fn}".R1.fastq.gz"
   ln -s ${fpath}".R2.fastq.gz" "/sc/arion/projects/psychAD/fastqs/"${2}"/Sample_"${l}/${fn}".R2.fastq.gz"
}
done
