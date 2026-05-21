#!/usr/bin/env python3

from __future__ import annotations
import argparse, subprocess
from pathlib import Path


def find_vcfs(input_dir: Path):
    return sorted(
        list(input_dir.glob("*.vcf")) +
        list(input_dir.glob("*.vcf.gz"))
    )


def build_lsf_script(input_dir: str, output_dir: str, plink_temp: str):

    return f"""#!/bin/bash

#BSUB -P acc_CommonMind
#BSUB -q express
#BSUB -o logs/%J_%I.out
#BSUB -e logs/%J_%I.err
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=3000]"
#BSUB -W 02:00

set -euo pipefail

INPUT_DIR="{input_dir}"
OUTPUT_DIR="{output_dir}"
PLINK_TEMP="{plink_temp}"

mkdir -p "$OUTPUT_DIR" logs "$PLINK_TEMP"

mapfile -t FILES < <(
  find "$INPUT_DIR" -maxdepth 1 -type f \\( -name "*.vcf" -o -name "*.vcf.gz" \\) | sort
)

NUM_FILES="${{#FILES[@]}}"

if [[ "$LSB_JOBINDEX" -gt "$NUM_FILES" ]]; then
  echo "[$(date)] Index exceeds file count"
  exit 0
fi

INPUT_FILE="${{FILES[$LSB_JOBINDEX-1]}}"

TEMP="${{INPUT_FILE##*/}}"
BASENAME="${{TEMP%.vcf.gz}}"
BASENAME="${{BASENAME%.vcf}}"

# ---- runtime mapping log ----
echo -e "${{LSB_JOBID}}\\t${{LSB_JOBINDEX}}\\t${{INPUT_FILE}}" >> logs/job_file_map.tsv

echo "Processing: $INPUT_FILE"

ml plink/1.90

plink --vcf "${{INPUT_FILE}}" \\
    --const-fid 0 \\
    --geno 0.10 \\
    --maf 0.10 \\
    --recode vcf \\
    --seed 42 \\
    --out "${{PLINK_TEMP}}/${{BASENAME}}"

ml bcftools/1.21

bcftools view -Oz \\
    --threads 12 \\
    -o "${{OUTPUT_DIR}}/${{BASENAME}}_UnFixedName.vcf.gz" \\
    "${{PLINK_TEMP}}/${{BASENAME}}.vcf"

# ---- fix sample name issue due to plink < v2 adding 0_ ----
bcftools query \\
    -l "${{OUTPUT_DIR}}/${{BASENAME}}_UnFixedName.vcf.gz" \\
    | awk 'BEGIN {{OFS="\t"}}{{c1=$1; sub(/^0_/, "", c1); print $1, c1}}' \\
    > "${{PLINK_TEMP}}/${{BASENAME}}_sampleNames.txt"

bcftools reheader \\
    -s "${{PLINK_TEMP}}/${{BASENAME}}_sampleNames.txt" \\
    -o "${{OUTPUT_DIR}}/${{BASENAME}}.vcf.gz" \\
    "${{OUTPUT_DIR}}/${{BASENAME}}_UnFixedName.vcf.gz"

bcftools index -t "${{OUTPUT_DIR}}/${{BASENAME}}.vcf.gz"

rm "${{OUTPUT_DIR}}/${{BASENAME}}_UnFixedName.vcf.gz"

"""


def submit_job(script_path: Path, n: int):
    cmd = [
        "bsub",
        "-J", f"filtVcf[1-{n}]",
        "<",
        str(script_path)
    ]

    # IMPORTANT: shell=True required for redirection syntax
    cmd_str = " ".join(cmd)

    print("\n Submitting job:")
    print(cmd_str)

    result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print("\n Submission failed:")
        print(result.stderr)
        raise SystemExit(1)

    print("\n Submission successful:")
    print(result.stdout.strip())


def get_argument_parser():
    """Generate and return argument parser."""

    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("-p", "--plink_temp", default="plink_temp")
    parser.add_argument("-q", "--queue", default="express")
    parser.add_argument("-n", "--dry-run", action="store_true")

    return parser


def main():
    """Main entry point"""

    parser = get_argument_parser()
    args = parser.parse_args()


    files = find_vcfs(args.input_dir)
    n = len(files)

    if n == 0:
        raise SystemExit("No VCF files found")

    script = build_lsf_script(args.input_dir, args.output_dir, args.plink_temp)

    script_path = Path("filter_perPool_vcf.bash")
    script_path.write_text(script)

    print(f"\n Script written: {script_path}")
    print(f" Files found: {n}")

    print("\n Submission command:")
    print(
        f"bsub -J 'filtVcf[1-{n}]' < {script_path}"
    )

    print("\n Runtime log will be written to:")
    print("logs/job_file_map.tsv")

    if args.dry_run:
        print("\n DRY RUN mapping:")
        for i, f in enumerate(files, 1):
            print(f"{i:4d} → {f}")
    else:
        print("\n Submitting the array job:")
        submit_job(script_path, n)


if __name__ == "__main__":
    main()
