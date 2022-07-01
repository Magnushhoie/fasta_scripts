# fasta_scripts
Collection of FASTA processing scripts

Usage:
```bash

# Split FASTA file into 3 entries of length 0-800, 800-1200 and 1200+ residues
fasta_split_fasta.py --fasta data/sequences.fasta

# Split FASTA file to single entry FASTA files
fasta_split_single_entries.py --fasta data/sequences.fasta --outdir output/

# Count duplicate sequences in FASTA
fasta_remove_duplicate_seqs.py --fasta data/sequences.fasta --output unique_seqs.fasta

```

Requirements:
```bash
# Install biopython
pip install -r requirements.txt
```
