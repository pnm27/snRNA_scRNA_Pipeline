#!/bin/bash


# ${1} is the donor name
# ${2} is the hash_file with donor name and the cell barcodes (with header, tab separated)
# ${3} is temp dir (with sample name as the parent dir)
# ${4} is the output dir (with sample name as the parent dir)
# ${5} is the multiplexed bam file

ml samtools

args=$@
mito_pref="MT"
verbose="FALSE"
POSITIONAL_ARGS=()
parse_args ()
{
    while [[ $# -gt 0 ]]
    do
    key="$1"

    case $key in
        -d|--donor)
        donor=$2
        shift
        shift
        ;;

        --hash_file)
        hash_file=$2
        shift
        shift
        ;;

        -t|--temp_dir)
        tempdir=$2
        shift
        shift
        ;;

        -o|--out_dir)
        outdir=$2
        shift
        shift
        ;;

        -b|--bam)
        pooled_bam=$2
        shift
        shift
        ;;

        -m|--mito_file)
        mito_file=$2
        shift
        shift
        ;;

        --bed)
        bedfile=$2
        shift
        shift
        ;;

        --mito_pref)
        mito_pref=$2
        shift
        shift
        ;;

        -v|--verbose)
        verbose="TRUE"
        shift
        ;;

        -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;

        *)
        POSITIONAL_ARGS+=("$1") # save positional arg
        shift
        ;;

    esac
    done
}

parse_args $args

# For backwards compatibility
set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters


if [[ $# -gt 0 ]]; then
    if [[ ! -z $1 ]]; then
        donor=${1}
    fi
    if [[ ! -z $2 ]]; then
        hash_file=${2}
    fi
    if [[ ! -z $3 ]]; then
        tempdir=${3}
    fi
    if [[ ! -z $4 ]]; then
        outdir=${4}
    fi
    if [[ ! -z $5 ]]; then
        pooled_bam=${5}
    fi
    if [[ ! -z $6 ]]; then
        mito_file=${6}
    fi
    if [[ ! -z $7 ]]; then
        bedfile=${7}
    fi
fi

# For some reason, awk wouln't work
if [ ! -d "${tempdir}" ]; then mkdir -p ${tempdir}; fi
if [ ! -d "${outdir}" ]; then mkdir -p ${outdir}; fi

if [[ ${verbose} == "TRUE" ]]; then
    set -x
fi

awk -v don="${donor}" '(NR> 1 && $1 == don){print $2}' ${hash_file} > ${tempdir}${donor}.txt
# samtools view -D CB:${3}${1}.txt ${5} -bho "${4}${1}.bam"
samtools view -D CB:${tempdir}${donor}.txt ${pooled_bam} -bho "${tempdir}${donor}.bam"
sleep 20
samtools index "${tempdir}${donor}.bam" &> /dev/null
sleep 20

if [[ ! -z ${mito_file} ]]; then
    mito_reads=$(samtools view -c "${tempdir}${donor}.bam" "${mito_pref}")
    echo "Number of mito reads for the donor ${donor}: ${mito_reads}" >> ${mito_file}
    # Remove unwanted contigs and secondary alignments
    samtools view -L "${bedfile}" -F 0x100 -o "${outdir}${donor}.bam" "${tempdir}${donor}.bam"
    sleep 20
    samtools index "${outdir}${donor}.bam" &> /dev/null && \
        ( rm "${tempdir}${donor}.bam" "${tempdir}${donor}.bam.bai"; \
        echo -n "Number of reads after filtering mito reads for the donor ${donor}: " >> ${mito_file}; \
            samtools view -c "${outdir}${donor}.bam" >> ${mito_file} ) && exit 0 || exit 1

else
    # Remove unwanted contigs and secondary alignments
    samtools view -L "${bedfile}" -F 0x100 -o "${outdir}${donor}.bam" "${tempdir}${donor}.bam"
    sleep 20
    samtools index "${outdir}${donor}.bam" &> /dev/null && \
        ( rm "${tempdir}${donor}.bam" "${tempdir}${donor}.bam.bai"; \
            samtools view -c "${outdir}${donor}.bam" >> ${mito_file} ) && exit 0 || exit 1

fi

if [[ ${verbose} == "TRUE" ]]; then
    set +x
fi
