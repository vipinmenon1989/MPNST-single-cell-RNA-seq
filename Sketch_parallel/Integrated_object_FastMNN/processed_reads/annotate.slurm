#!/bin/bash
#SBATCH --job-name=S_1M_Analysis
#SBATCH --account=ihc
#SBATCH --partition=ihc        # Use a high-memory node
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=16         # Parallelize clustering
#SBATCH --mem=800G                 # 500GB is safer for 1M cell SNN graphs
#SBATCH --time=224:00:00
#SBATCH --output=annotate_1M_%j.log

cd /local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/Integrated_object_FastMNN/processed_reads/

source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

Rscript heatmaps.R
