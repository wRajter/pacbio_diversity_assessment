#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=6
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=log_output_%A.log

# Creating a reference tree based on the reference alignment using RAxML

# Variables
REF_ALIGNMENT="reference_data/reference_alignment_gblocksStrict.fasta"
OTU_DIR="phyl_placement/Gblocks_strict"
MODEL="GTR+G"
NCORES=6
PREFIX='T2_strict'

# Activate conda environment
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate phylo_placment

mkdir -p ${OTU_DIR}/

# Compute tree
raxml-ng \
    --msa ${REF_ALIGNMENT} \
    --model ${MODEL} \
    --prefix ${PREFIX} \
    --threads ${NCORES}


mv ${PREFIX}* ${OTU_DIR}/
