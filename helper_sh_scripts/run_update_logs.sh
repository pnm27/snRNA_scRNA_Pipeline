#!/usr/bin/bash
# Runs helper_py_scripts/update_logs.py for the samples in INPUT_FILE
# Load necessary Modules

declare -A params # contains the parameter for python script
declare -A args # contains the values for above params

# DON'T CHANGE THE BELOW SET
params["BAM_DIR"]="-b"
params["PICARD_DIR"]="-p"
params["MAP_FILE"]="-m"
params["OUTPUT_FILE"]="-o"
params["BAM_STRUCT"]="--bam_struct"
params["PICARD_STRUCT"]="--pc_struct"
params["DEMUX_STRUCT"]="--dem_struct"
params["DEMUX_SUFFIX"]="--dem_info"
params["WET_LAB_FILE"]="-w"
params["ANNOT"]="--common_annotations"
params["SWAP_CORRECT_DF"]="--swap_correct"
params["DONOR_NAME_CONV_DF"]="--donorName_correct"

# Outputs that update_logs.py can compile (see its --help)
declare -A opt_flags=(
    [ss_l]="STARsolo Log.final.out"
    [ss_g_f]="STARsolo Gene Features.stats"
    [ss_gf_f]="STARsolo GeneFull Features.stats"
    [ss_g_s]="STARsolo Gene Summary.csv"
    [ss_gf_s]="STARsolo GeneFull Summary.csv"
    [ss_bc]="STARsolo Barcodes.stats"
    [pc_gc]="Picard GC summary metrics"
    [pc_rs]="Picard RNAseq metrics"
)

REQUIRED_ARGS=("BAM_DIR" "OUTPUT_FILE" "BAM_STRUCT")


# CHANGE THIS ACCORDING TO PROJECT REQUIREMENTS
INPUT_FILE="/sc/arion/projects/psychAD/pnm/fastq_files.txt" # REQUIRED
args["BAM_DIR"]="/sc/arion/projects/psychAD/STARsolo_bams/" # REQUIRED
args["PICARD_DIR"]="/sc/arion/projects/psychAD/STARsolo_bams/"
args["MAP_FILE"]="/sc/arion/projects/psychAD/pnm/Final_out_MAP_2.tsv"
args["OUTPUT_FILE"]="/sc/arion/projects/psychAD/pnm/All_logs.tsv" # REQUIRED
args["BAM_STRUCT"]="Sample_<sample>*/" # REQUIRED
args["PICARD_STRUCT"]="Sample_<sample>*/"
args["DEMUX_STRUCT"]="Sample_<sample>*/"
args["DEMUX_SUFFIX"]="_STARsolo_vS_info.tsv"
args["WET_LAB_FILE"]=""
args["ANNOT"]=""
args["SWAP_CORRECT_DF"]=""
args["DONOR_NAME_CONV_DF"]=""

# Outputs to compile and their file suffixes. Leave a value empty to use
# the default suffix of update_logs.py; remove (or comment out) a line to
# leave that output out.
declare -A opt_values=(
    [ss_l]=""
    [ss_g_f]=""
    [ss_gf_f]=""
    [ss_g_s]=""
    [ss_gf_s]=""
    [ss_bc]=""
    [pc_gc]=""
    [pc_rs]=""
)
verbosity="" # Can be set to -v, -vv or -vvv for verbose output, else leave it empty


# Fails fast if a required arg (or the input file) is missing
check_required_args() {
    local missing=()
    local key

    for key in "${REQUIRED_ARGS[@]}"; do
        [[ -z "${args[$key]}" ]] && missing+=("${key}")
    done

    if (( ${#missing[@]} > 0 )); then
        echo "ERROR: required args not set: ${missing[*]}" >&2
        exit 1
    fi

    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "ERROR: INPUT_FILE '${INPUT_FILE}' doesn't exist" >&2
        exit 1
    fi
}


# Adds --<key> [suffix] for each output in opt_values; unknown keys are
# skipped with a warning
append_optional_flags() {
    local key

    for key in "${!opt_values[@]}"; do
        if [[ -z "${opt_flags[$key]+set}" ]]; then
            echo "WARNING: '${key}' is not a known output of update_logs.py, skipping it" >&2
            continue
        fi
        cmd_args+=("--${key}")
        [[ -n "${opt_values[$key]}" ]] && cmd_args+=("${opt_values[$key]}")
    done
}


check_required_args

samp_list=()
while IFS=$'\t' read -r -a l
do
{
    if [[ -n "${l[0]}" ]] && ! grep -q "^#" <<< "${l[0]}"; then
        # Remove 'cDNA' at the end of the names
        # glob can add cDNA or HTO as required
        samp_name=$(rev <<< "${l[0]}" | cut -d "/" -f1 | rev | sed "s#[-_][cC][dD][nN][aA]##g" )
        samp_list+=("${samp_name}")
    fi
}
done < "${INPUT_FILE}"

if (( ${#samp_list[@]} == 0 )); then
    echo "ERROR: no samples found in INPUT_FILE '${INPUT_FILE}'" >&2
    exit 1
fi

# build command
cmd_args=()
for key in "${!args[@]}"
do
    [[ -n "${args[$key]}" ]] && cmd_args+=("${params[$key]}" "${args[$key]}")
done
append_optional_flags
[[ -n "${verbosity}" ]] && cmd_args+=("${verbosity}")

echo "STARTING the script ${0} at:"
date
echo
set -x
# Restrict glob expansion
python3 "$(dirname "$0")/../helper_py_scripts/update_logs.py" "${cmd_args[@]}" \
    --samples "${samp_list[@]}"


set +x
echo
echo "FINISHED the script at:"
date
