#!/bin/bash
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-h200-1
#SBATCH --job-name=seurat_fastMNN_1M
#SBATCH --output=logs/fastMNN_mnn_%j.out
#SBATCH --error=logs/fastMNN_mnn_%j.err
#SBATCH --time=96:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=450G  

# -------------------------------------------------------------------------
# LOGGING AND ENVIRONMENT SETUP
# -------------------------------------------------------------------------
echo "--- JOB STARTED AT $(date) ---"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"

# Load Conda Environment
# Ensure this environment has 'batchelor' and 'SeuratWrappers' installed
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

# -------------------------------------------------------------------------
# EXECUTION
# -------------------------------------------------------------------------
# Using the single-script approach for better memory management
Rscript FastMN.R

echo "--- JOB FINISHED AT $(date) ---"
