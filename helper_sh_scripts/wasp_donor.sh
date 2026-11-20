#!/usr/bin/env bash
# The steps of a donor bam's WASP job that are the same for every modality
# (rules/wasp/rna.smk, rules/wasp/atac.smk); only the remapping between
# them differs (STAR, bwa-mem). As github.com/bmvdgeijn/WASP's mapping
# pipeline; its scripts are in $WASP_HOME (the WASP env module), run with
# that module's python.
#
#   prepare  the bam and the genotypes have chromosomes in common; the
#            reads over the donor's SNPs (to remap, as <prefix>.remap.fq*.gz)
#            and the others (<prefix>.keep.bam), with <prefix> =
#            <work>/<bam name without .bam>
#   finish   keeps the remapped reads that map back to the same place,
#            merged with the others, sorted and indexed; optionally without
#            the reads marked as duplicates (0x400); flagstats and the
#            donor's stats row; removes <work>
#
# Usage:
#   wasp_donor.sh prepare --work DIR --bam BAM --samples "G1 G2" \
#       --snp_index H5 --snp_tab H5 --haps H5 [--paired] --log LOG
#   wasp_donor.sh finish --work DIR --bam BAM --remapped BAM --out BAM \
#       --threads N --filter_log LOG --intersect_log LOG \
#       --flagstat_before F --flagstat_after F [--flagstat_wasp F] \
#       [--drop_duplicates] [--paired] --pool P --donor D --final_donor F \
#       --sample "G1;G2" \
#       --modality cDNA --stats TSV

set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${WASP_HOME:?WASP_HOME is not set: load the WASP module (references: tools: wasp)}"

cmd=${1:?prepare or finish}
shift
paired=""
drop_duplicates=""
flagstat_wasp=""
threads=1
while [[ $# -gt 0 ]]; do
    case $1 in
        --paired) paired=1; shift ;;
        --drop_duplicates) drop_duplicates=1; shift ;;
        --work) work=$2; shift 2 ;;
        --bam) bam=$2; shift 2 ;;
        --samples) samples=$2; shift 2 ;;
        --snp_index) snp_index=$2; shift 2 ;;
        --snp_tab) snp_tab=$2; shift 2 ;;
        --haps) haps=$2; shift 2 ;;
        --log) log=$2; shift 2 ;;
        --remapped) remapped=$2; shift 2 ;;
        --out) out=$2; shift 2 ;;
        --threads) threads=$2; shift 2 ;;
        --filter_log) filter_log=$2; shift 2 ;;
        --intersect_log) intersect_log=$2; shift 2 ;;
        --flagstat_before) flagstat_before=$2; shift 2 ;;
        --flagstat_after) flagstat_after=$2; shift 2 ;;
        --flagstat_wasp) flagstat_wasp=$2; shift 2 ;;
        --pool) pool=$2; shift 2 ;;
        --donor) donor=$2; shift 2 ;;
        --final_donor) final_donor=$2; shift 2 ;;
        --sample) sample=$2; shift 2 ;;
        --modality) modality=$2; shift 2 ;;
        --stats) stats=$2; shift 2 ;;
        *) echo "wasp_donor.sh: unknown option $1" >&2; exit 1 ;;
    esac
done
prefix="${work}/$(basename "${bam}" .bam)"

case ${cmd} in
prepare)
    rm -rf "${work}"
    mkdir -p "${work}"
    # One genotype sample ID per line (find_intersecting_snps.py's --samples)
    # shellcheck disable=SC2086
    printf '%s\n' ${samples} > "${work}/samples.txt"

    # The bam and the genotypes have chromosomes in common (as WASP matches
    # them: as they are, or with the 'chr' prefix added or stripped)
    samtools view -H "${bam}" \
        | awk -F'\t' '$1 == "@SQ" {for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) print substr($i, 4)}' \
        > "${work}/contigs.txt"
    python "${REPO}/helper_py_scripts/wasp_check_inputs.py" contigs \
        --contigs "${work}/contigs.txt" --snp_tab "${snp_tab}" --bam "${bam}"

    # Reads over the donor's SNPs (to remap) and the others (kept)
    python "${WASP_HOME}/mapping/find_intersecting_snps.py" ${paired:+--is_paired_end} \
        --is_sorted --snp_index "${snp_index}" --snp_tab "${snp_tab}" \
        --haplotype "${haps}" --samples "${work}/samples.txt" \
        --output_dir "${work}" "${bam}" &> "${log}"
    ;;

finish)
    # Keep the reads that map back to the same place, with those over no SNP
    python "${WASP_HOME}/mapping/filter_remapped_reads.py" "${prefix}.to.remap.bam" \
        "${remapped}" "${work}/remap.keep.bam" &> "${filter_log}"
    samtools merge -f -@ "${threads}" "${work}/merged.bam" \
        "${prefix}.keep.bam" "${work}/remap.keep.bam"
    if [[ -n ${drop_duplicates} ]]; then
        # Duplicates as marked by the aligner (cellranger-arc: per cell barcode)
        samtools sort -@ "${threads}" -T "${work}/sort" -o "${work}/sorted.bam" "${work}/merged.bam"
        samtools view -@ "${threads}" -b -F 0x400 -o "${out}" "${work}/sorted.bam"
    else
        samtools sort -@ "${threads}" -T "${work}/sort" -o "${out}" "${work}/merged.bam"
    fi
    samtools index -@ "${threads}" "${out}"

    samtools flagstat -@ "${threads}" "${bam}" > "${flagstat_before}"
    samtools flagstat -@ "${threads}" "${out}" > "${flagstat_after}"
    wasp_arg=()  # (expanded so that an empty array is fine with set -u in bash < 4.4)
    if [[ -n ${flagstat_wasp} ]]; then
        if [[ -n ${drop_duplicates} ]]; then
            samtools flagstat -@ "${threads}" "${work}/sorted.bam" > "${flagstat_wasp}"
        else
            cp "${flagstat_after}" "${flagstat_wasp}"
        fi
        wasp_arg=(--flagstat_wasp "${flagstat_wasp}")
    fi
    python "${REPO}/helper_py_scripts/wasp_stats.py" parse --pool "${pool}" \
        --donor "${donor}" --final_donor "${final_donor:-${donor}}" --modality "${modality}" --sample "${sample}" \
        ${paired:+--paired_end} --intersect_log "${intersect_log}" \
        --filter_log "${filter_log}" --flagstat_before "${flagstat_before}" \
        --flagstat_after "${flagstat_after}" ${wasp_arg[@]+"${wasp_arg[@]}"} -o "${stats}"
    rm -rf "${work}"
    ;;

*)
    echo "wasp_donor.sh: unknown step ${cmd} (prepare or finish)" >&2
    exit 1
    ;;
esac
