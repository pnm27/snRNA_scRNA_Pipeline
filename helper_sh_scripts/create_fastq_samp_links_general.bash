#!/usr/bin/bash

###
#Create input files called r1_samples, r2_samples, etc
# cd fastq_dir/
# find -name "*R1*fastq.gz" | sort > projdir/all_r1_fastqs
# cd projdir
# sed -i 's#^\.#projdir#' all_r1_fastqs

#Run:
# sh create_fastq_samp_links_general.sh all_r1_fastqs 

# Creating input for Snakemake's run
# If R1 files have paths like:
# fastq_dir/myFile_A_L001_R1_001.fastq.gz
# we want to retain only myFile-A
# rev all_r1_fastqs  | cut -d "/" -f1 | cut -d "_" -f4- | rev | sed "s#_#-#" | sort > fastq_files.txt

echoinfo ()
{
    echo -e "\033[32m[INFO] $@\033[0m"
}

help_msg ()
{
    echoinfo "This script is used to create fastq_lengths fastq files. "
}

usage_msg ()
{
    echoinfo "create_fastq_samp_links_general.bash -d <output_directory> -i <r1_fastq_file> [-n] [-h]"
    echoinfo ""
    echoinfo "-d | --link_dir     Provide an output directory where the fastq links will be created. Default ./fastq_links"
    echoinfo "-i | --input          File containing absolute paths to 'R1' file for each sample/pool. REQUIRED"
    echoinfo "-n | --dry_run        A dry-run which will run a max of 5 lines from the input file."
    echoinfo ""
    echoinfo "NOTE: While using dry-run, although ONLY the first 5 lines will be read the output may be on more "
    echoinfo "than 1 line as the script considers each line as a sample and thus, the I1, I2, I3, R1 and R2 files "
    echoinfo "will be searched and executed accordingly."
}

args=$@
link_dir="./fastq_links"
dry_run="false"
input_file=""
parse_args ()
{
    while [[ $# -gt 0 ]]
    do
        key="$1"
        case $key in
            -d|--link_dir)
            link_dir="$2"
            shift
            shift
            ;;

            -i|--input)
            input_file="$2"
            shift
            shift
            ;;

            -n|--dry_run)
            dry_run="true"
            shift
            ;;

            -h|--help)
            help_msg
            usage_msg
            exit 0
            ;;

            -*|--*)
            echo "Unknown option $1"
            exit 1
            ;;

        esac
    done
}

parse_args $args

###
# HARCODED: DEPRACATED STYLE
# link_dir="fastq_links"
# If relative path is given then change it to an absolute path
link_dir=$(realpath "${link_dir}" )

(
if [[ ${dry_run} == "true" ]]
then
    head -5 "${input_file}"
else
    cat "${input_file}"
fi
) | while read l
do
{
   # BD2-Set-12-cDNA-a_S1_L001_R1_001.fastq.gz
   # BD2-Set-5-ATAC-b_S1_L001_R1_001.fastq.gz
   # base_filename=BD2-Set-12-a, BD2-Set-5-b
   # Expected structure
   # ATAC/file_dir/files
   # cDNA/file_dir/files

   # Project-based base_filename selection
   # base_filename=$(rev <<< ${l} | cut -d "/" -f1 | cut -d "_" -f5- | rev | sed "s#[aA][tT][aA][cC]-##g" | sed "s#[cC][dD][nN][aA]-##g")
   # base_filename=$(sed "s#NW_#NW-#g" <<< ${base_filename} | sed "s#NY_#NY-#g" | sed "s#_#-#g")
   base_filename=$(rev <<< ${l} | cut -d "/" -f1 | cut -d "_" -f4- | rev | sed "s#[aA][tT][aA][cC]-##g" | sed "s#[cC][dD][nN][aA]-##g" | sed "s#_#-#g" | sed -- "s#-cDNA##g" )

   # subdir=$(echo ${l} | tr '[:upper:]' '[:lower:]' | grep -q "atac" && echo "ATAC" || echo "cDNA") # OLD STYLE
   subdir=$(echo ${l} | tr '[:upper:]' '[:lower:]') # TEMP Store
   if grep -q -E "atac" <<< ${subdir}; then
       subdir="ATAC"
   elif grep -q -E "hto" <<< ${subdir}; then
       subdir="HTO"
   else
       subdir="cDNA"
   fi
   # subdir="cDNA"
   # echo "var value: ${subdir}"
   ext=""
   # If the filename ends in "gz" pick the next term for extension otherwise just the last string following
   # "."
   r3_file=""
   i2_file=""
   if [[ ${l##*gz} == "" ]]; then
      ext=$(rev <<< ${l} | cut -d "." -f1-2 | rev)
   else
      ext=$(rev <<< ${l} | cut -d "." -f1 | rev)
   fi
   # Extract Lane numbers (for samples sequenced across multiple lanes)
   # lane_num=$(sed -r "s#.*_(L[0-9]+)_.*#\1#g" <<< "${l}")
   r1_suffix=$(rev <<< ${l} | cut -d "/" -f1 | cut -d "_" -f-3 | rev )
   r1_suffix=${r1_suffix%${ext}}
   # r1_suffix="R1"
   r2_file=$(sed -r "s#([._]*)R1([._]*)#\1R2\2#g" <<< "${l}")
   r2_suffix=$(sed -r "s#([._]*)R1([._]*)#\1R2\2#g" <<< "${r1_suffix}")
   i1_file=$(sed -r "s#([._]*)R1([._]*)#\1I1\2#g" <<< "${l}")
   i1_suffix=$(sed -r "s#([._]*)R1([._]*)#\1I1\2#g" <<< "${r1_suffix}")
   if [[ ${subdir} == "ATAC" ]]; then
      r3_file=$(sed -r "s#([._]*)R1([._]*)#\1R3\2#g" <<< "${l}")
      r3_suffix=$(sed -r "s#([._]*)R1([._]*)#\1R3\2#g" <<< "${r1_suffix}")
      [[ ${dry_run} == "true" ]] && echo "Create dir: ${link_dir%/}/${subdir}/${base_filename}" || mkdir -p "${link_dir%/}/${subdir}/${base_filename}"
   else
      i2_file=$(sed -r "s#([._]*)R1([._]*)#\1I2\2#g" <<< "${l}")
      i2_suffix=$(sed -r "s#([._]*)R1([._]*)#\1I2\2#g" <<< "${r1_suffix}")
      [[ ${dry_run} == "true" ]] && echo "Create dir: ${link_dir%/}/${subdir}/${base_filename}" || mkdir -p "${link_dir%/}/${subdir}/${base_filename}"
   fi

   
   # Link R1 file
   [[ ${dry_run} == "true" ]] && echo "R1: ${l} ${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${r1_suffix}${ext}" || ln -s "${l}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${r1_suffix}${ext}"

   # Check if R2 file exists
   if [ -f "${r2_file}" ]; then    
      [[ ${dry_run} == "true" ]] && echo "R2: ${r2_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${r2_suffix}${ext}" || ln -s "${r2_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${r2_suffix}${ext}"
   else
      echo "Couldn't find R2 file; Linked the R1 file only!"
   fi

   # Check if I1 file exists
   if [ -f "${i1_file}" ]; then    
      [[ ${dry_run} == "true" ]] && echo "I1: ${i1_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${i1_suffix}${ext}" || ln -s "${i1_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${i1_suffix}${ext}"
   else
      echo "Couldn't find I1 file!"
   fi

   if [[ ${subdir} == "ATAC" ]]; then
      # Check if R3 file exists
      if [ -f "${r3_file}" ]; then    
         [[ ${dry_run} == "true" ]] && echo "R3: ${r3_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${r3_suffix}${ext}" || ln -s "${r3_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${r3_suffix}${ext}"
      else
         echo "Couldn't find ${r3_file} file!"
      fi
   else
      # Check if I2 file exists
      if [ -f "${i2_file}" ]; then    
         [[ ${dry_run} == "true" ]] && echo "I2: ${i2_file}" "${link_dir%/}/${base_filename}/${base_filename}_${i2_suffix}${ext}" || ln -s "${i2_file}" "${link_dir%/}/${subdir}/${base_filename}/${base_filename}_${i2_suffix}${ext}"
      else
         echo "Couldn't find ${i2_file} file!"
      fi
   fi

}
done
