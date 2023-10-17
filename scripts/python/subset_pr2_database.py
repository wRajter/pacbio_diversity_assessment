from Bio import SeqIO
import re
import os

def subset_sequences_by_taxon(input_file_path, output_file_path, tax_level='f', limit=1, max_length=3000):
    """
    Subset sequences from a given input fasta file based on a specified taxonomic level, sequence count limit,
    and maximum sequence length.

    Parameters:
    - input_file_path (str): Path to the input fasta file.
    - output_file_path (str): Path to the output fasta file.
    - tax_level (str): Taxonomic level (e.g., 'f' for family, 'g' for genus). Default is 'f'.
    - limit (int): Maximum sequences for one taxonomic level. Default is 1.
    - max_length (int): Maximum sequence length. Default is 3000.

    Returns:
    - None: Writes the subsetted sequences to the output fasta file.
    """

    # Initialize dictionary to hold sequences per taxon level
    tax_dict = {}

    # Read and parse sequences
    with open(input_file_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # Check the length of the sequence
            if len(record.seq) <= max_length:
                # Extract taxonomic information from header
                tax_match = re.search(rf'{tax_level}:([^,]+)', record.description)
                if tax_match:
                    taxon = tax_match.group(1)

                    # Add sequence to tax_dict, if taxon count is below the threshold/limit
                    if taxon not in tax_dict:
                        tax_dict[taxon] = [record]
                    elif len(tax_dict[taxon]) < limit:
                        tax_dict[taxon].append(record)

    # Write filtered sequences to a new file
    with open(output_file_path, 'w') as output_handle:
        for taxon_records in tax_dict.values():
            SeqIO.write(taxon_records, output_handle, "fasta")
