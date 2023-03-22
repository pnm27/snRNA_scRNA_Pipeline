#!/usr/bin/bash
# spreadsheet_dir is the directory with all "xlsx" file present under each preparer's dir i.e. Aram will have her own folder, so will Clara, Stathis, etc.
# spreadsheet_file is the file that will have all data compiled for use by snakemake pipeline
# converter_file was provided by Jaro and contains 'set' vs 'subID' values
# run_script is the script that will be executed based on certain conditions (explained just before it)
# process_log will record which files has been processed into the "spreadsheet file"
# files_tracker_log will contain information related to each file and the samples present in it related to rerunning required or not
# i.e. if a sample's or multiple samples in a file has been updated wrt its name or HTO numbers or the barcodes associated with them
# It needs to be removed from the "spreadsheet file" and be updated with the newer one. It will also serve the purpose to instruct the user
# to rerun those samples for downstream processing

spreadsheet_dir="/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/"
spreadsheet_file="/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/New_wet_lab_spreadsheet.tsv"
converter_file="/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/converter.xlsx"
run_script="/sc/arion/projects/psychAD/pnm/helper_py_scripts/create_wet_lab_info.py"
process_log="/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/process_log.txt"
files_tracker_log="/sc/arion/projects/psychAD/pnm/NPSAD_spreadsheet/files_tracker_log.txt"
all_f="" # To store file name(s)
proj="NPSAD" # This name will be searched in wet lab files' column
# columns="unique_sample_ID hashtag ab_barcode SubID Set_number" # Columns that are relevant to downstream analyses

echo "Starting the script $0 at: "
date

echo "Creating process log (if it doesn't exist already, just 'touch'ing the log file): "
if [ ! -s "${process_log}" ]; then touch ${process_log}; fi

declare -i c=0

# If the files are present inside subfolders
# Inside the dir (mentioned in spreadsheet_dir)
for var in `ls ${spreadsheet_dir%/}/*/*.xlsx`
do
    # File is not empty (doesn't matter if it exists but empty) or file name not present in "proccess log" file
    # or xlsx file is newer than the spreadsheet
    if [ ! -s "${spreadsheet_file}" ] || ! grep -q "${var}" "${process_log}" || [ "${var}" -nt "${spreadsheet_file}" ]
    then
        c+=1
        all_f="${all_f} ${var}"

    elif grep -q "${var}" "${process_log}"
    then
        echo "The file: ${var} has been already processed"
    elif [ "${var}" -ot "${spreadsheet_file}" ]
    then
        echo "The file, ${var}, is older than the output file-${spreadsheet_file}"
    else
        echo "Irrelevant condition"
    fi
done
# Remove leading and trailing spaces, merge continuous spaces into one
all_f=$(echo "${all_f}" | xargs)

if [ ! -z "${all_f}" ]
then
    if [ ! -z "${converter_file}" ]
    then
        echo "Processing ${c} files, with a converter file, starting at: "
        date
        python3 ${run_script} ${all_f} -o ${spreadsheet_file} -c ${converter_file} -l ${files_tracker_log} -p ${proj} # --columns ${columns}
    else
        echo "Processing ${c} files, without a converter file, starting at: "
        date
        python3 ${run_script} ${all_f} -o ${spreadsheet_file} -l ${files_tracker_log} -p ${proj}
    fi
fi

read -a l <<< "${all_f}"
for i in "${l[@]}"
do
    echo "${i}" >> ${process_log}
done

echo "Finished processing ${c} files at: "
date