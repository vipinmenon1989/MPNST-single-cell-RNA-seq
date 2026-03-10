# MPNST-single-cell-RNA-seq

Single-cell RNA sequencing (scRNA-seq) analysis pipeline for **Malignant Peripheral Nerve Sheath Tumor (MPNST)** — a rare and aggressive soft-tissue sarcoma. This repository contains R scripts for quality control, normalization, batch integration, and dimensionality reduction of multi-sample MPNST scRNA-seq datasets.

## Overview

MPNST is a highly malignant tumor arising from peripheral nerve sheaths, commonly associated with Neurofibromatosis type 1 (NF1). This pipeline processes aggregated single-cell transcriptomic data from multiple batches/labs, performing rigorous QC filtering and RPCA-based batch integration to enable unbiased cell-type discovery and comparative analysis across samples.

The pipeline is structured as a sequential two-step workflow:

1. **Part 1 — QC & Preprocessing** (`MPNST.R`): Load CellBender-filtered data, apply quality control, normalize, find variable features, and run PCA
2. **Part 2 — Batch Integration** (`MPNST_RPCA_Integration.R`): Perform RPCA (Reciprocal PCA) integration across batches and generate UMAP visualizations

Additional sketch-based approximate analysis scripts are provided for rapid exploratory analysis of large datasets.

## Repository Structure

```
MPNST-single-cell-RNA-seq/
├── MPNST.R                    # Part 1: QC, normalization, PCA (outputs seurat_obj_pre_integration.rds)
├── MPNST.slurm                # SLURM job script for Part 1
├── MPNST_RPCA_Integration.R   # Part 2: RPCA batch integration + UMAP (30-40 hr runtime)
├── RPCA.slurm                 # SLURM job script for Part 2
├── Sketch.R                   # Sketch-based approximate analysis (single core)
├── Sketch.slurm               # SLURM job for Sketch.R
├── Sketch_parallel.R          # Parallelized sketch analysis
├── Sketch_parallel/           # Output directory for sketch results
```

## Pipeline Details

### Part 1: QC and Preprocessing (`MPNST.R`)

**Input**: `aggregated_cellbender_filtered.rds` — a merged Seurat object with CellBender-corrected counts from multiple samples

**Steps**:
1. HPC parallel setup (4 workers, 200 GB memory limit)
2. Pre-filter QC scatter plots (UMI counts vs mitochondrial %, counts vs genes)
3. Pre-filter summary statistics (grouped by batch and lab)
4. Quality filtering: remove cells with mitochondrial fraction > 20%
5. Post-filter QC visualization
6. Normalization with `NormalizeData`
7. Highly variable feature selection (`FindVariableFeatures`, 5000 features, VST method)
8. Data scaling (`ScaleData`)
9. PCA (30 principal components)

**Output**: `seurat_obj_pre_integration.rds` — checkpoint file for Part 2

### Part 2: RPCA Batch Integration (`MPNST_RPCA_Integration.R`)

**Input**: `seurat_obj_pre_integration.rds`

**Steps**:
1. Sequential mode (1 core, ~125 GB RAM — avoids memory spikes)
2. UMAP on unintegrated data (for before/after comparison)
3. RPCA integration via Seurat's `IntegrateLayers` (k.weight = 100)
4. UMAP on integrated data

**Estimated Runtime**: 30–40 hours on HPC

**Output**: Integrated Seurat object with UMAP embeddings and QC plots

## Usage

### Running on HPC Cluster (SLURM)

```bash
# Step 1: QC and preprocessing
sbatch MPNST.slurm

# Step 2: Batch integration (run after Step 1 completes)
sbatch RPCA.slurm
```

### Running Locally

```r
# Part 1
source("MPNST.R")

# Part 2 (after Part 1 completes)
source("MPNST_RPCA_Integration.R")
```

### Sketch-Based Approximate Analysis (optional, for rapid exploration)

```bash
sbatch Sketch.slurm
```

## Input Data Requirements

- **Format**: Seurat RDS object (`aggregated_cellbender_filtered.rds`)
- **Source**: CellBender-filtered count matrices aggregated across samples
- **Metadata**: Must include `orig.ident` (batch/sample column) and `lab` columns

## Output Files

| File | Description |
|------|-------------|
| `seurat_obj_pre_integration.rds` | Post-QC Seurat object ready for integration |
| `01_QC_Pre_Scatter_Uniform.png` | Pre-filter QC scatter plots |
| `02_QC_Post_Scatter_Uniform.png` | Post-filter QC scatter plots |
| `03_UMAP_Unintegrated.png` | UMAP of unintegrated data |
| `04_UMAP_Integrated.png` | UMAP after RPCA integration |
| `pre_filter_summary.csv` | Pre-filter cell count statistics per batch |
| `post_filter_summary.csv` | Post-filter cell count statistics per batch |

## R Dependencies

```r
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr"))
install.packages("future")  # For parallel processing
```

Seurat v5 is required (uses the split-layer architecture for batch integration).

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `nfeatures` | 5000 | Number of highly variable features |
| `npcs` | 30 | PCA dimensions |
| `percent.mt` cutoff | 20% | Mitochondrial fraction filter |
| `k.weight` | 100 | RPCA integration anchor weight |
| `dims` | 1:30 | Dimensions used for UMAP |
| `workers` | 4 (Part 1), 1 (Part 2) | Parallel workers |

## Author

**Vipin Menon**
PhD Candidate, Computational Biology (Genome Editing)
Hanyang University, Seoul, South Korea
GitHub: [vipinmenon1989](https://github.com/vipinmenon1989)
