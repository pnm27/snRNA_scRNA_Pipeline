#!/bin/bash
#BSUB -P acc_CommonMind
#BSUB -q express
#BSUB -o logs/%J_%I.out
#BSUB -e logs/%J_%I.err
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=3000]"
#BSUB -W 01:00

set -euo pipefail

# -------- HELP FUNCTION --------
usage() {
  cat <<'EOF'
Usage: filter_perPool_vcf.bash <INPUT_DIR> <OUTPUT_DIR> [plink_temp]

Description:
  Run an LSF array job where each task processes one VCF file
  from INPUT_DIR and writes results to OUTPUT_DIR.

Arguments:
  INPUT_DIR   Directory containing .vcf files
  OUTPUT_DIR  Directory where output will be written
  plink_temp  Optional dir to store temp plink outputs

Example:
  N=\$(find data/vcfs -maxdepth 1 -name "*.vcf.gz" | wc -l)
  bsub -J "filtVcf[1-\${N}]" filter_perPool_vcf.bash data/vcfs results/output

Notes:
  - Make sure N matches the number of input files
  - Logs are written to ./logs/
EOF
}

# -------- ARGUMENT PARSING --------
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$#" -lt 2 ]]; then
  echo "Error: Missing required arguments"
  usage
  exit 1
fi


INPUT_DIR="" # argparse required argument
OUTPUT_DIR="" # argparse required argument
PLINK_TEMP="plink_temp" # argparse optional argument

PLINK_TEMP="${PLINK_TEMP%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"

# -------- VALIDATION --------
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: INPUT_DIR does not exist: $INPUT_DIR"
  exit 1
fi

mkdir -p "$OUTPUT_DIR" logs "${PLINK_TEMP}"

# -------- BUILD FILE LIST --------
mapfile -t FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | sort)

NUM_FILES="${#FILES[@]}"

# -------- ARRAY OVERSHOOT PROTECTION --------
if [[ "$LSB_JOBINDEX" -gt "$NUM_FILES" ]]; then
  echo "[$(date)] Task $LSB_JOBINDEX exceeds file count ($NUM_FILES). Exiting."
  exit 0
fi

INPUT_FILE="${FILES[$LSB_JOBINDEX-1]}"
# BASENAME=$(basename "$INPUT_FILE" .vcf.gz)
# Remove path
TEMP="${INPUT_FILE##*/}"
# Remove .vcf.gz or .vcf
BASENAME="${TEMP%.vcf.gz}"
BASENAME="${BASENAME%.vcf}"

echo "[$(date)] Task ${LSB_JOBINDEX} processing: ${INPUT_FILE}"

ml plink/1.90

# -------- MAIN COMMAND --------
plink --vcf "${INPUT_FILE}" \
    --geno 0.10 \
    --maf 0.10 \
    --recode vcf \
    --seed 42 \
    --out "{PLINK_TEMP}/${BASENAME}"


ml bcftools/1.21

bcftools view -W'tbi' \
    --threads=12 \
    -Oz \
    -o "${OUTPUT_DIR}/${BASENAME}.vcf.gz" \
    "{PLINK_TEMP}/${BASENAME}.vcf"