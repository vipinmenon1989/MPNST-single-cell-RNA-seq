#!/bin/bash
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-h200-1
#SBATCH --job-name=seurat_fastMNN_1M_cr
#SBATCH --output=/local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/cr_filtered/logs/cr_fastMNN_mnn_%j.out
#SBATCH --error=/local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/cr_filtered/logs/cr_fastMNN_mnn_%j.err
#SBATCH --time=200:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=450G  

# -------------------------------------------------------------------------
# LOGGING AND ENVIRONMENT SETUP
# -------------------------------------------------------------------------
echo "--- JOB STARTED AT $(date) ---"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"

cd /local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/cr_filtered/
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
