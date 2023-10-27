#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=6
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=T2_Moderate_log_output_%A.log

# Creating a reference tree based on the reference alignment using RAxML

# Variables
REF_ALIGNMENT="pacbio_diversity_assessment/alignments/pr2_subset/reference_alignment_gblocksModerate.fasta"
OTU_DIR="pacbio_diversity_assessment/reference_trees/pr2_subset/gblocksModerate"
MODEL="GTR+FO+G4m"
NCORES=6
PREFIX='T2_Moderate'

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