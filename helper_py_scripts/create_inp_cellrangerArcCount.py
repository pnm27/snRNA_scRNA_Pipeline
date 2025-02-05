#!/usr/bin/env python3

"""Creates a list of effective barcodes for cellSNP run.

This script takes a count matrix as an input and creates a list of 
'effective' barcodes, dependent on the condition if some cells should 
be removed or not (e.g. if already a demultiplexing tool was run with an
h5ad column containing classifications and if some of the cells that
have been classified as a doublet and/or negative cells have to be removed 
then too this script can be used).

"""


import pandas as pd
import os, argparse
# using datetime module


def get_argument_parser():
    """Generate and return argument parser."""

    # Parse Command-Line arguments
    parser = argparse.ArgumentParser(description="Produce metadata "
    "input for cellranger arc count. Takes all the fastq files as input "
    "and classifies them as cDNA or ATAC by infering the path."
    )

    parser.add_argument('fastqs', nargs='+', help="Path to the fastq "
    "files, including both modals - cDNA and ATAC.",
    )

    # Optional parameters
    parser.add_argument('-o', '--output', 
    help="Name of the output file",
    default="metadata.csv"
    )
    parser.add_argument('-s', '--sample_id', 
    help="Sample ID value",
    default="sample", dest='sample'
    )
    
    return parser


def main():
    """Main function"""

    parser = get_argument_parser()
    args = parser.parse_args()
    cdna_files = []
    atac_files = []
    for fastq in args.fastqs:
        if 'atac' in fastq.lower():
            atac_files.append(fastq)
        elif 'cdna' in fastq.lower():
            cdna_files.append(fastq)
        else:
            raise SyntaxError("Couldn't figure out the modality of the "
                f"file from the filename: {fastq}")
        
    cdna_pat = [ os.path.dirname(f) for f in cdna_files ]
    atac_pat = [ os.path.dirname(f) for f in atac_files ]

    # Sanity check
    assert len(set(cdna_pat)) == 1, "All cDNA files don't share a parent dir!"
    assert len(set(atac_pat)) == 1, "All ATAC files don't share a parent dir!"

    op_dict = {
        "fastqs": [cdna_pat[0], atac_pat[0]],
        "sample": [args.sample, args.sample],
        "library_type": ["Gene Expression", "Chromatin Accessibility"]
    }
    op_df = pd.DataFrame(op_dict)
    op=args.output if args.output.endswith('.csv') else args.output + '.csv'
    op_df.to_csv(op, index=False, sep=',')


if __name__ == '__main__':
    main()