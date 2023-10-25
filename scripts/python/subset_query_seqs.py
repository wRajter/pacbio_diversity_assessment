# Imports
import os
from Bio import SeqIO


# Functions

# Extract sequence IDs with Bacteria or Archaea in their taxonomic assignment
def find_bacteria_archaea(tax_assign_dir, file):
    '''
    Extract sequence IDs with Bacteria or Archaea in their taxonomic assignment.

    Parameters:
    - tax_assign_dir (str): Directory path containing taxonomic assignment files.
    - file (str): Taxonomic assignment file.

    Returns:
    - bacterial_and_archaeal_ids (set): Set of sequence IDs with either Bacteria or Archaea in their taxonomic assignment.
    '''
    bacterial_and_archaeal_ids = []
    with open(os.path.join(tax_assign_dir, file), "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            taxonomy = parts[1]
            seq_id = parts[0]
            if 'k:Bacteria' in taxonomy or 'k:Archaea' in taxonomy:
                bacterial_and_archaeal_ids.append(seq_id)
    return bacterial_and_archaeal_ids

# Filter fasta files
def filter_fasta_file(input_file, output_file, ids_to_exclude):
    '''
    Filters sequences in a fasta file based on a list of sequence IDs to exclude.

    Parameters:
    - input_file (str): Path to the input fasta file.
    - output_file (str): Path where the filtered fasta file should be saved.
    - ids_to_exclude (set): Set of sequence IDs to exclude.
    '''
    sequences = [record for record in SeqIO.parse(input_file, "fasta") if record.id not in ids_to_exclude]
    SeqIO.write(sequences, output_file, "fasta")
