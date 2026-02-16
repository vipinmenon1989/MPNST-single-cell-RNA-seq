#!/bin/bash
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-h200-1
#SBATCH --job-name=seurat_rpca
#SBATCH --output=logs/seurat_rpca_sketch_parallel_%j.out
#SBATCH --error=logs/seurat_rpca_sketch_parallel_%j.err
#SBATCH --time=96:00:00 # Reduced time; this will now finish much faster
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=450G  # Increased to ensure the 430GB object fits comfortably

# -------------------------

echo "--- JOB STARTED AT $(date) ---"

# Load Conda Environment
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "Running on node: $(hostname)"
echo "Requested CPUs: $SLURM_CPUS_PER_TASK"
echo "Requested Memory: $SLURM_MEM_PER_NODE"

# EXECUTE THE MASTER SCRIPT
# Ensure 'main_pipeline.R' and 'integrate_sketch.R' are in this folder.
Rscript main_pipeline.R

echo "--- JOB FINISHED AT $(date) ---"