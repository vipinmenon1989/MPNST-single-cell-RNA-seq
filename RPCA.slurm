#!/bin/bash
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --job-name=seurat_rpca
#SBATCH --output=logs/seurat_rpca_%j.out    # Standard output log
#SBATCH --error=logs/seurat_rpca_%j.err     # Error log
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G

# -------------------------

echo "--- JOB STARTED AT $(date) ---"
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "Running on node: $(hostname)"
echo "Requested CPUs: $SLURM_CPUS_PER_TASK"
echo "Requested Memory: $SLURM_MEM_PER_NODE"

Rscript MPNST_RPCA_Integration.R

echo "--- JOB FINISHED AT $(date) ---"