# split_fasta.py
# Magnus Haraldson Høie, maghoi.dtu.dk June 2022

import argparse
import logging
import os
from collections import defaultdict

from Bio import SeqIO

usage = """split_fasta.py -f long_fasta.fasta --limits '0,800,1200'
Outputs 3 FASTA files with sequence lengths 0-800, 800-1200, and 1200+"""
parser = argparse.ArgumentParser(
    description="Split FASTA file by sequence lengths", usage=usage
)
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA to split")
parser.add_argument(
    "-l", "--limits", default="0,800,1200", help="FASTA length limits, comma-separated"
)
parser.add_argument("-o", "--output", default="output/", help="Output directory")
dict_args = vars(parser.parse_args())

try:
    dict_args["limits"] = dict_args["limits"].split(",")
    dict_args["limits"] = list(map(int, dict_args["limits"]))
except:
    print(
        f"ERROR: option --limits must be comma separated positive integers (e.g. 0,800,1200)"
    )

print(dict_args)


def get_fasta_stats(fasta_dict):
    """Calculates distribution of sequences in FASTA dict"""

    sequences = [v.seq for v in fasta_dict.values()]
    lengths = [len(s) for s in sequences]

    log.info(f"{len(fasta_dict)} IDs, {len(sequences)} sequences")
    log.info(f"Max length {max(lengths)}, min length {min(lengths)}")


def parse_fasta_file(fasta_path, verbose=False):
    """Reads in FASTA file, prints statistics"""

    fasta_path = "data/sequences.fasta"
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    log.info(f"Read {len(fasta_dict)} FASTA entries from {fasta_path}")
    get_fasta_stats(fasta_dict)

    return fasta_dict


def check_overlapping_sequences(fasta_dict, verbose=False):
    """Checks which IDs are overlapping in FASTA dict"""

    inverseDict = defaultdict(list)

    for k, v in fasta_dict.items():
        seq = str(v.seq)
        inverseDict[seq].append(k)

    for k in list(inverseDict.keys()):
        if len(inverseDict[k]) < 2:
            del inverseDict[k]

    log.info(f"Found overlapping sequences for {len(inverseDict)} IDs")

    if verbose:
        for k, v in inverseDict.items():
            log.info(f"{v} : {k}")

    return inverseDict


def filter_fasta_dict_on_seqlen(fasta_dict, min_len=0, max_len=800):
    """Filters FASTA dict based on sequence length"""

    out_dict = {
        k: v
        for k, v in fasta_dict.items()
        if len(v.seq) >= min_len and len(v.seq) < max_len
    }
    log.info(
        f"Filtering to {len(out_dict)} / {len(fasta_dict)} sequences >= {min_len} and < {max_len}"
    )

    return out_dict


def save_fasta_dict(outfile, fasta_dict):
    """Saves FASTA dict to file"""

    if len(fasta_dict) >= 1:
        log.info(f"Saving to {outfile}\n")
        with open(outfile, "w") as out_handle:
            SeqIO.write(fasta_dict.values(), out_handle, "fasta")


def main(dict_args):
    """Split FASTA"""

    fasta_path = dict_args["fasta"]
    os.makedirs("output", exist_ok=True)

    # Read in FASTA to dict
    fasta_dict = parse_fasta_file(fasta_path)

    # Check for overlap
    _ = check_overlapping_sequences(fasta_dict)

    # Save filtered FASTA files
    limits = dict_args["limits"]
    for i, min_len in enumerate(limits):
        if i + 1 <= len(limits) - 1:
            max_len = limits[i + 1]
        else:
            max_len = int(10e4)

        # Write file
        filtered_dict = filter_fasta_dict_on_seqlen(
            fasta_dict, min_len=min_len, max_len=max_len
        )
        save_fasta_dict(f"output/lengths_{min_len}-{max_len}.fasta", filtered_dict)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[{asctime}] {message}", style="{")
    log = logging.getLogger(__name__)
    log.info("Splitting FASTA file")

    main(dict_args)
