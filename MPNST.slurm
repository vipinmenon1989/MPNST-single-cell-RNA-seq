#!/bin/bash
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-h200-1
#SBATCH --job-name=seurat_rpca
#SBATCH --output=logs/seurat_%j.out    # Standard output log
#SBATCH --error=logs/seurat_%j.err     # Error log
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=300G

# -------------------------

echo "--- JOB STARTED AT $(date) ---"
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "Running on node: $(hostname)"
echo "Requested CPUs: $SLURM_CPUS_PER_TASK"
echo "Requested Memory: $SLURM_MEM_PER_NODE"

Rscript MPNST.R

echo "--- JOB FINISHED AT $(date) ---"