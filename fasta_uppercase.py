#!/Users/maghoi/opt/anaconda3/envs/py39/bin/python
# fasta_uppercase.py: Converts FASTA sequences to uppercase

import argparse
import logging
from pathlib import Path

# Fasta entries
from Bio import SeqIO


def cmdline_args():
    # Make parser object
    p = argparse.ArgumentParser(
        description="""
    Converts entries in input FASTA file to individual output FASTA files
    """
    )

    p.add_argument(
        "--fasta",
        required=True,
        help="Input fasta file to split into individual fasta entries",
    )
    p.add_argument(
        "--outfile", required=False, default="", help="Output directory"
    )

    return p.parse_args()


def main(fasta_path, outfile):
    """
    Splits entries in FASTA file into individual entry FASTA output files
    NB: Converts sequences to upper-case
    """

    log.info(f"Saving to {outfile}")

    # Read in FASTA entries
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    # Upper-case
    for k in fasta_dict.keys():
        fasta_dict[k].seq = fasta_dict[k].seq.upper()

    # Write as uppercase
    with open(outfile, "w") as out_handle:
        SeqIO.write(fasta_dict.values(), out_handle, "fasta")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[{asctime}] {message}", style="{")
    log = logging.getLogger(__name__)

    args = cmdline_args()
    log.info(f"Converting {args.fasta} to uppercase in {args.outfile} ...")

    # Add _upper to outfile unless specified
    if len(args.outfile) == 0:
        stem = Path(args.fasta).stem
        outfile = f"{stem}_upper.fasta"

    # Main
    main(args.fasta, outfile)
