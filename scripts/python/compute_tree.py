# Imports
import subprocess
import os


def run_raxmlng_check(ref_alignment, output_dir, model='GTR+G', prefix="T1", log_file='log_file.txt'):
    '''
    Run raxml-ng to check the alignment and save the output in the specified directory.
    '''
    # Specify command for raxm-ng
    cmd = [
        'raxml-ng', '--check',
        '--msa', ref_alignment,
        '--model', model,
        '--prefix', prefix
    ]
    # Run the command via the subrocess module
    results = subprocess.run(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                check=True)

    # Print output and error messages (if any)
    print(results.stdout)
    if results.stderr:
        print(results.stderr)


    # Move files to the output directory
    for file in os.listdir():
        if file.startswith(prefix):
            os.rename(file, os.path.join(output_dir, file))
    if results.returncode == 0:
        print(f'Output files were moved to: {output_dir}')

    return results.returncode



def generate_slurm_script(inputs):
    '''
    Generate the slurm script based on user inputs.
    '''
    script_content = f"""#!/bin/bash
#SBATCH --mail-user={inputs['email']}
#SBATCH --mail-type=END
#SBATCH --cpus-per-task={inputs['ncores']}
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output={inputs['prefix']}_log_output_%A.log

# Creating a reference tree based on the reference alignment using RAxML

# Variables
REF_ALIGNMENT="{inputs['ref_alignment']}"
OTU_DIR="{inputs['otu_dir']}"
MODEL="{inputs['model']}"
NCORES={inputs['ncores']}
PREFIX='{inputs['prefix']}'

# Activate conda environment
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate phylo_placment

mkdir -p ${{OTU_DIR}}/

# Compute tree
raxml-ng \\
    --msa ${{REF_ALIGNMENT}} \\
    --model ${{MODEL}} \\
    --prefix ${{PREFIX}} \\
    --threads ${{NCORES}}

mv ${{PREFIX}}* ${{OTU_DIR}}/"""

    return script_content



def submit_sbatch(ssh_client, slurm_script_path, output_file_path):
    """
    Submit a slurm script via SSH and save the resulting job ID to a file.

    Args:
    - ssh_client (paramiko.SSHClient): Active SSH client.
    - slurm_script_path (str): Path to the slurm script on the remote server.
    - output_file_path (str): Path to save the job ID on the local machine.

    Returns:
    - str: Job ID if submission was successful, None otherwise.
    """

    # Execute sbatch command
    command = f'sbatch {slurm_script_path}'
    stdin, stdout, stderr = ssh_client.exec_command(command)

    # Grab the output and error
    output = stdout.read().decode().strip()
    error = stderr.read().decode().strip()

    # Print only if there's meaningful content
    if output:
        print("Output:", output)
    if error:
        print("Error:", error)

    # Extracting the job ID from the output
    job_id = None
    if "Submitted batch job" in output:
        job_id = output.split()[-1]

    # Saving the job ID to the provided path
    if job_id:
        with open(output_file_path, 'w') as file:
            file.write(job_id)
    else:
        print("Failed to retrieve the job ID!")
        return None

    return job_id
