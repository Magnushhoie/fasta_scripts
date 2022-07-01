#!/Users/maghoi/opt/anaconda3/envs/py39/bin/python
# split_fasta.py: Converts entries in input FASTA file to individual output FASTA files
# By Magnus H. Høie, maghoi@dtu.dk, March 2022

import argparse
import logging
import os
import re
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
    p.add_argument("--outdir", required=True, default="", help="Output directory")

    return p.parse_args()


def split_fasta(fasta_path, outdir_path):
    """
    Splits entries in FASTA file into individual entry FASTA output files
    NB: Converts sequences to upper-case
    """

    # Make outdir
    os.makedirs(outdir_path, exist_ok=True)
    log.info(f"Saving to {outdir_path}")

    # Read in FASTA entries
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    # Write individual sequence FASTA files
    for id, seqrec in fasta_dict.items():
        id_new = re.search("\w+", id).group(0)

        if id_new is None:
            log.error(f"Error: Broken id {id_new} from {id} for sequence {seqrec}")
            raise Exception

        outfile = Path(outdir_path, id_new + ".fasta")

        if os.path.exists(outfile):
            log.error(f"Error: {outfile} already exists!")
        else:
            with open(outfile, "w") as out_handle:
                # Upper case conversion to not break colabfold_batch
                seqrec.seq = seqrec.seq.upper()
                SeqIO.write(seqrec, out_handle, "fasta")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[{asctime}] {message}", style="{")
    log = logging.getLogger(__name__)

    args = cmdline_args()
    log.info(f"Splitting entries in {args.fasta}, saving to {args.outdir} ...")

    split_fasta(args.fasta, args.outdir)
