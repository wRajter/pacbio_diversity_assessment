# Imports
from Bio import SeqIO
import re
import logging
import os
import subprocess
from Bio import AlignIO
import pandas as pd
import shutil

# Global variables
GBLOCKS_PATH = os.path.join('..', 'raw_data', 'packages', 'Gblocks_0.91b', 'Gblocks')


def subset_sequences_by_taxon(input_file_path, output_file_path,
                              tax_level='f', limit=1, min_length= 1500, max_length=2000,
                              log_file='subset_sequences.log'):
    '''
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
    '''

    # Configure logging
    log_path = os.path.join(os.path.dirname(output_file_path), log_file)
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[logging.FileHandler(log_path), logging.StreamHandler()])

    logging.info('Starting subsetting sequences...')
    logging.info(f'Using input file: {input_file_path}')
    logging.info(f'Writing to output file: {output_file_path}')
    logging.info(f'Taxonomic level: {tax_level}')
    logging.info(f'Limit per taxonomic group: {limit}')
    logging.info(f'Maximum sequence length: {max_length}')

    # Initialize dictionary to hold sequences per taxon level
    tax_dict = {}
    total_sequences = 0
    retained_sequences = 0

    # Read and parse sequences
    with open(input_file_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            total_sequences += 1
            # Logging every 1000 sequences for progress tracking
            if total_sequences % 1000 == 0:
                logging.info(f'Processed {total_sequences} sequences...')
            # Check the length of the sequence
            if min_length <= len(record.seq) <= max_length:
                # Extract taxonomic information from header
                tax_match = re.search(rf'{tax_level}:([^,]+)', record.description)
                if tax_match:
                    taxon = tax_match.group(1)

                    # Add sequence to tax_dict, if taxon count is below the threshold/limit
                    if taxon not in tax_dict:
                        tax_dict[taxon] = [record]
                        retained_sequences += 1
                    elif len(tax_dict[taxon]) < limit:
                        tax_dict[taxon].append(record)
                        retained_sequences += 1

    logging.info(f'Total sequences processed: {total_sequences}')
    logging.info(f'Sequences retained: {retained_sequences}')

    # Write filtered sequences to a new file
    with open(output_file_path, 'w') as output_handle:
        for taxon_records in tax_dict.values():
            SeqIO.write(taxon_records, output_handle, 'fasta')

    logging.info(f'Finished subsetting sequences. Check {output_file_path} for results.')




def run_mafft(input_file, output_file, log_file):
    '''Run MAFFT alignment on the input file using the FFT-NS-2 strategy.'''

    cmd = ['mafft', '--retree', '2', input_file]

    with open(output_file, 'w') as out_f, open(log_file, 'w') as log_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=log_f, text=True, check=True)

    if result.returncode == 0:
        print(f'Alignment was saved to: {output_file}')
        print(f'Log file was saved to: {log_file}')

    return result.returncode

def alignment_stats(alignment_file):
    # Load the alignment
    alignment = AlignIO.read(alignment_file, "fasta")

    # Number of sequences
    num_sequences = len(alignment)

    # Alignment length
    alignment_length = alignment.get_alignment_length()

    if alignment_length == 0:
        raise ValueError(f"The alignment in {alignment_file} has a length of zero. Please check the file.")

    # Calculate the total number of gaps across all sequences
    total_gaps = sum(str(record.seq).count('-') for record in alignment)

    # Average number of gaps per sequence
    average_gaps = total_gaps / num_sequences

    # Percentage of gaps across the sequences
    percentage_gaps = (average_gaps / alignment_length) * 100

    # Percentage identity (for simplicity, we'll just calculate this for pairwise comparisons)
    total_columns = alignment_length * num_sequences
    identical_columns = sum(col.count(col[0]) for col in zip(*alignment))
    percentage_identity = (identical_columns / total_columns) * 100

    return num_sequences, alignment_length, average_gaps, percentage_gaps, percentage_identity





def creating_df_from_alignment(alignment_path):
    '''
    Extract taxonomy data from sequence headers in a FASTA alignment and create a DataFrame.

    This function assumes that the sequence headers in the FASTA file are in the format:
    >[sequence_ID]_U;tax=k:[kingdom],d:[domain],p:[phylum],c:[class],o:[order],f:[family],g:[genus],s:[species]

    Parameters:
    - alignment_path (str): Path to the FASTA alignment file.

    Returns:
    - DataFrame: A pandas DataFrame containing the taxonomy data extracted from the FASTA file. Columns are:
        - ID: Sequence identifier
        - kingdom: Kingdom level taxonomy
        - domain: Domain level taxonomy
        - phylum: Phylum level taxonomy
        - class: Class level taxonomy
        - order: Order level taxonomy
        - family: Family level taxonomy
        - genus: Genus level taxonomy
        - species: Species level taxonomy
    '''

    taxopaths = []

    # Open the input file for reading and output file for writing
    with open(alignment_path, 'r') as infile:
        # Read through each line in the input file
        for line in infile:
            # If the line starts with '>', it's a header line
            if line.startswith('>'):
                taxopaths.append(line)


    # List to hold dictionaries for DataFrame creation
    data = []

    # Process each taxopath
    for taxopath in taxopaths:
        # Splitting at underscore to separate ID and taxonomic path
        ID, taxonomy = taxopath.split(';tax=')
        # Removing '>'
        ID = ID[1:]

        # Dictionary to store taxonomic data for this sequence ID
        taxo_data = {'ID': ID}

        # Split taxonomy string into taxonomic ranks and extract relevant information
        for rank_info in taxonomy.split(','):
            # Using partition to split at the first occurrence of ':'
            rank, _, name = rank_info.partition(':')
            taxo_data[rank] = name.strip()

        # Append the dictionary to the list
        data.append(taxo_data)

    # Convert list of dictionaries into a DataFrame
    df = pd.DataFrame(data)

    # Rearrange columns in the desired order
    cols_order = ['ID', 'k', 'd', 'p', 'c', 'o', 'f', 'g', 's']
    df = df[cols_order]

    # Rename columns
    df.columns = ['ID', 'kingdom', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    return df

def subset_fasta_based_on_dataframe(fasta_path, output_path, dataframe):
    '''Subset a fasta file based on a Pandas dataframe.'''
    ids_to_keep = dataframe['ID'].tolist()
    ids_set = set(ids_to_keep)
    write_seq = False
    with open(fasta_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                seq_id = line.split(';tax=')[0][1:]  # Extract sequence ID from header
                if seq_id in ids_set:
                    write_seq = True
                    outfile.write(line)
                else:
                    write_seq = False
            elif write_seq:
                outfile.write(line)

def run_gblocks(input_alignment, b1, b2, b3, b4, b5, name_specifier='', gblocks_path = GBLOCKS_PATH):

    # Variables
    directory_path = os.path.dirname(input_alignment)
    input_file_name, _ = os.path.splitext(os.path.basename(input_alignment))
    # Specify the original and desired output names
    original_fasta_output = input_alignment + ".gb"
    fasta_output = os.path.join(directory_path, f'{input_file_name}_gblocks{name_specifier}.fasta')
    original_html_output = input_alignment + ".gb.htm"
    html_output = os.path.join(directory_path, f'{input_file_name}_gblocks{name_specifier}.fasta.html')
    # Path for the output .log file
    log_file_path = os.path.join(directory_path, f'{input_file_name}_gblocks{name_specifier}.log')

    # Gblocks parameters
    cmd = [gblocks_path,
           input_alignment,
           '-t=d',
           f'-b1={b1}',
           f'-b2={b2}',
           f'-b3={b3}',
           f'-b4={b4}',
           f'-b5={b5}',
           '-e=.gb']

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Rename the files
    os.rename(original_fasta_output, fasta_output)
    os.rename(original_html_output, html_output)

    # Grab info from the gblocks HTML file and save it to the log file
    # Open the HTML file and read its content
    with open(html_output, 'r') as file:
        content = file.read()

        # Using regex to find content between <pre> and </pre>
        match = re.search(r'<pre>.*?<b>Parameters used</b>(.*?)</pre>', content, re.DOTALL)

        # Start with the Gblocks stdout content
        log_content = f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}\n\n\nParameters used:\n"

        if match:
            extracted_info = match.group(1).strip()  # Extract the content
            log_content += extracted_info  # Append the extracted info to the log content
        else:
            log_content += "Pattern not found in the HTML file!"

        # Write the combined content to the .log file
        with open(log_file_path, 'w') as log_file:
            log_file.write(log_content)


def run_gblocks_grid_search(mafft_alignment,
                            b1_values,
                            b2_values,
                            b3_values,
                            b4_values,
                            b5_values,
                            gblocks_path):
    '''
    Execute Gblocks using a grid of parameter settings to fine-tune the alignment process.

    Parameters:
    - mafft_alignment (str): Path to the alignment file generated using MAFFT.
    - b1_values (list): List of values for b1 parameter.
    - b2_values (list): List of values for b2 parameter.
    - b3_values (list): List of values for b3 parameter.
    - b4_values (list): List of values for b4 parameter.
    - b5_values (list): List of values for b5 parameter.
    - gblocks_path (str): Path to the Gblocks executable.

    Returns:
    None. However, the function will produce alignment files with various parameters settings based on the grid search.
    '''
    for b5_value in b5_values:
        for i in range(len(b1_values)):
            b1 = b1_values[i]
            b2 = b2_values[i]
            b3 = b3_values[i]
            b4 = b4_values[i]
            b5 = b5_value
            run_gblocks(input_alignment = mafft_alignment,
                        b1 = b1,
                        b2 = b2,
                        b3 = b3,
                        b4 = b4,
                        b5 = b5,
                        name_specifier=f'_grid_{b5_value}_{i}',
                        gblocks_path = gblocks_path)



def move_search_grid_files(directory_path, grid_search_path):
    '''
    Transfer Gblocks output files to a designated 'grid_search' directory and delete all associated .html files.

    Parameters:
    - directory_path (str): Path to the directory containing the output files from Gblocks.
    - grid_search_path (str): Destination path for the 'grid_search' directory where output files will be moved.

    Returns:
    None. However, the function will move the specified files and delete any .html files in the target directory.
    '''
    # Move all the output files to the separate 'grid_search' directory
    files_to_move = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if 'grid' in f.lower() and os.path.isfile(os.path.join(directory_path, f))]

    grid_search_path = os.path.join(directory_path, 'grid_search')

    # Check if directory already exists, if not, create it
    if not os.path.exists(grid_search_path):
        os.mkdir(grid_search_path)

    for file in files_to_move:
        file_name = os.path.basename(file)  # Extract just the filename from the full path
        destination_path = os.path.join(grid_search_path, file_name)  # Construct the destination path
        shutil.move(file, destination_path)

    # delete all the html files
    for file in os.listdir(grid_search_path):
        if file.endswith('.html'):
            file_path = os.path.join(grid_search_path, file)
            os.remove(file_path)
