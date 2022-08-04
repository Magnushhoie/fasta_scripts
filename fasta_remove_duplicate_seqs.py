#!/home/maghoi/.conda/envs/py37/bin/python
# Removes duplicate sequences in FASTA file
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Removes duplicate seqs from input file, writes to input_out")
parser.add_argument("-f", "--fasta", dest="INPUT_FASTA", required=True, help="Input fasta file")
parser.add_argument("-o", "--ouput", dest="OUTPUT_FASTA", help="Output fasta file")
args = parser.parse_args()  

if args.OUTPUT_FASTA:
    output_fasta = args.OUTPUT_FASTA
else:
    output_fasta = args.INPUT_FASTA + ".filtered"

with open(output_fasta, 'w') as outFile:
    dict_seqs = {}
    dict_ids = {}
    count_entries = 0

    for record in SeqIO.parse(args.INPUT_FASTA, 'fasta'):
        count_entries += 1

        if str(record.seq) not in dict_seqs:
            # Store unique sequences
            dict_seqs[str(record.seq)] = True

            # Write ID+sequence if sequence is unique
            SeqIO.write(record, outFile, 'fasta')

            if str(record.id) not in dict_ids:
                # Store unique IDs
                dict_ids[str(record.id)] = True
            
# Print statistics
print("Total entries:", count_entries)
print("Unique IDs:", len(dict_ids.keys()))
print("Unique sequences:", len(dict_seqs.keys()))
print("Removed sequences in %s: %s" % (output_fasta, count_entries - len(dict_seqs.keys())))
