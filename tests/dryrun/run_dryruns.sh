#!/usr/bin/env bash
# Dry runs (snakemake -n) of the pipeline for each scenario in scenarios/.
# Builds a throwaway project (empty fastqs and reference files) from
# new_config.example.yaml, so no data or cluster is needed.
#
# Usage: bash tests/dryrun/run_dryruns.sh [scenario ...]
#   SNAKEMAKE  snakemake command (default: snakemake)
#   SNAKEFILE  Snakefile to run (default: <repo>/Snakefile)
#   KEEP_TMP   set to keep the throwaway project
#
# Each scenario is a config override with header comments:
#   # expect: <rules that must be in the DAG>
#   # absent: <rules that must not be in the DAG>
#   # resources: <rule> <key=value ...>   (resolved resources of that rule)
#   # path: <text>                        (must be in an input/output path)
#   # fails: <text>                       (snakemake must fail with this text)
#   # donors: <donor ...>                 (the jobs up to the donor checkpoint
#                                          'create_inp_splitBams' are done and
#                                          found these donors, so per-donor
#                                          jobs are in the DAG too)
#
# Each scenario runs in its own copy of the throwaway project.

set -uo pipefail

TEST_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "${TEST_DIR}/../.." && pwd)
SNAKEMAKE=${SNAKEMAKE:-snakemake}
SNAKEFILE=${SNAKEFILE:-${REPO}/Snakefile}

TMP=$(mktemp -d)
[[ -z "${KEEP_TMP:-}" ]] && trap 'rm -rf "${TMP}"' EXIT

# Throwaway project ---------------------------------------------------------
mkdir "${TMP}/base"
cp "${REPO}/new_config.example.yaml" "${TMP}/base/new_config.yaml"
cd "${TMP}/base"
mkdir -p projdir/fastq_files/Sample_P1-cDNA projdir/fastq_files/Sample_P1-HTO \
    prev/cellSNP/Sample_P1-cDNA
printf "P1-cDNA\n" > projdir/fastq_files.txt
for m in cDNA HTO; do
    for r in R1 R2; do
        touch "projdir/fastq_files/Sample_P1-${m}/P1-${m}_S1_L001_${r}_001.fastq.gz"
    done
done
printf "SetNumber\tGt_count\nP1\t4\n" > projdir/pool_info.tsv
printf "SetNumber\tGt_count\tdonors\tpool_vcf\nP1\t4\ta,b,c,d\tprojdir/P1.vcf.gz\n" \
    > projdir/pool_info4.tsv
touch projdir/P1.vcf.gz projdir/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.sorted.vcf.gz \
    prev/cellSNP/Sample_P1-cDNA/cellSNP.cells.vcf.gz
# Outputs of an earlier hashsolo run (run_tag wo_gt) and a swap correction file
mkdir -p projdir/final_count_matrix/solo/wo_gt/P1 projdir/demultiplex/info/solo/wo_gt/P1
touch projdir/final_count_matrix/solo/wo_gt/P1/P1-HTO_STARsolo_out.h5ad \
    projdir/demultiplex/info/solo/wo_gt/P1/P1-HTO_hashsolo_stats.tsv
printf "pool\tversion\tfrom\tto\nP1\tsolo:wo_gt\tdonor0\tD1\n" > projdir/swap.tsv
# identify_swaps
mkdir -p outdir && touch outdir/reference.vcf.gz


dry_run() {
    PYTHONPATH="${REPO}" ${SNAKEMAKE} -s "${SNAKEFILE}" -n \
        --configfile new_config.yaml "$1" 2>&1
}


# Creates the outputs of the jobs planned up to the donor checkpoint and
# writes its donor hash file with 250 cells for each donor
mark_checkpoint_done() {
    local sc=$1 f h d i
    shift
    dry_run "${sc}" | sed -n 's/^    output: //p' | tr ',' '\n' | sed 's/^ *//' \
        | grep -v '^resolved/' | while read -r f; do
            [[ -z "${f}" ]] && continue
            mkdir -p "$(dirname "${f}")" && touch "${f}"
        done
    sleep 1 # hash files newer than their inputs
    find . -name '*_bc_hash.txt' | while read -r h; do
        printf "Subj_ID\tbarcodes\n" > "${h}"
        for d in "$@"; do
            for i in $(seq 250); do printf "%s\tBC%s_%s\n" "${d}" "${d}" "${i}"; done
        done >> "${h}"
    done
}

# Scenarios -----------------------------------------------------------------
scenarios=("$@")
if (( ${#scenarios[@]} == 0 )); then
    scenarios=("${TEST_DIR}"/scenarios/*.yaml)
else
    scenarios=("${scenarios[@]/#/${TEST_DIR}/scenarios/}")
    scenarios=("${scenarios[@]/%.yaml/}")
    scenarios=("${scenarios[@]/%/.yaml}")
fi

n_fail=0
for sc in "${scenarios[@]}"; do
    name=$(basename "${sc}" .yaml)
    cp -r "${TMP}/base" "${TMP}/${name}"
    cd "${TMP}/${name}"

    donors=$(sed -n 's/^# donors://p' "${sc}")
    [[ -n "${donors}" ]] && mark_checkpoint_done "${sc}" ${donors}

    out=$(dry_run "${sc}")
    status=$?
    errors=()

    fails=$(sed -n 's/^# fails: *//p' "${sc}")
    if [[ -n "${fails}" ]]; then
        (( status == 0 )) && errors+=("snakemake didn't fail")
        grep -qF "${fails}" <<< "${out}" || errors+=("no '${fails}' in the output")
    elif (( status != 0 )); then
        errors+=("snakemake exited with ${status}")
    fi

    jobs=$(grep -oE "^(local)?(rule|checkpoint) [A-Za-z0-9_]+" <<< "${out}" \
        | awk '{print $2}' | sort -u)
    for r in $(sed -n 's/^# expect://p' "${sc}"); do
        grep -qx "${r}" <<< "${jobs}" || errors+=("missing rule ${r}")
    done
    for r in $(sed -n 's/^# absent://p' "${sc}"); do
        grep -qx "${r}" <<< "${jobs}" && errors+=("unexpected rule ${r}")
    done
    while read -r rule kvs; do
        [[ -z "${rule}" ]] && continue
        res=$(awk -v r="${rule}" '$0 ~ "^(local)?rule "r":" {f=1} f && /^    resources:/ {print; exit}' <<< "${out}")
        for kv in ${kvs}; do
            grep -qE "(^|[ ,])${kv}(,|$)" <<< "${res}" || errors+=("${rule}: expected ${kv} in '${res# *}'")
        done
    done < <(sed -n 's/^# resources://p' "${sc}")
    for p in $(sed -n 's/^# path://p' "${sc}"); do
        grep -qE "^    (input|output|log): .*${p}" <<< "${out}" || errors+=("no path with '${p}'")
    done

    if (( ${#errors[@]} == 0 )); then
        echo "PASS  ${name}"
    else
        n_fail=$(( n_fail + 1 ))
        echo "FAIL  ${name}"
        printf '      %s\n' "${errors[@]}"
        (( status != 0 )) && [[ -z "${fails}" ]] && tail -n 15 <<< "${out}" | sed 's/^/      | /'
    fi
done

echo
echo "${n_fail} of ${#scenarios[@]} scenarios failed"
(( n_fail == 0 ))
