# fasta_scripts
Collection of FASTA processing scripts

Usage:
```bash

# Split FASTA file into 3 entries of length 0-800, 800-1200 and 1200+ residues
src/split_fasta.py --fasta data/sequences.fasta

# Split FASTA file to single entry FASTA files
split_fasta_single_entries.py --fasta data/sequences.fasta --outdir output/

```

Requirements:
```bash
# Install biopython
pip install -r requirements.txt
```
