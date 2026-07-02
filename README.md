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
├── Sketch.R                   # Sketch-based approximate analysis (single core, in-process)
├── Sketch.slurm               # SLURM job for Sketch.R
├── Sketch_parallel.R          # Earlier in-process-parallel draft of Sketch.R (kept for reference)
├── Sketch_parallel/           # Extended FastMNN / sketch-integration pipeline (scripts + SLURM jobs,
│                               # NOT an output directory -- see "Extended Pipeline" below)
│   ├── main_pipeline.R        # Recommended sketch+RPCA integration (hands off to integrate_sketch.R)
│   ├── integrate_sketch.R     # Sub-process launched by main_pipeline.R for parallel RPCA
│   ├── FastMN.R                # Alternative integration via FastMNN (SeuratWrappers/batchelor)
│   ├── r_metadata.R / bp.R / bp_export.R / bp_cell_converter.R / R_escape.R
│   │                           # Export helpers (BPCells / metadata / Python-bridge formats)
│   ├── cr_filtered/            # Same FastMNN pipeline, alternate CellRanger-filtered input
│   └── Integrated_object_FastMNN/
│       ├── diag.R / unitegrated.R / Plots.R / Immune_plots.R / integrated_plot.R / split_lab.R
│       ├── lab_subsets/wu_lab.R
│       └── processed_reads/    # Downstream clustering, annotation, marker/heatmap generation
├── Snakefile                  # Snakemake workflow tying the above stages together
├── config/config.yaml         # Snakemake configuration (input paths, resources)
├── envs/environment.yaml      # Conda environment (R + Seurat + friends) for --use-conda
└── .github/workflows/ci.yml   # CI: R syntax checks + Snakemake dry-run
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

### Extended Pipeline: FastMNN Sketch Integration (`Sketch_parallel/`)

For very large datasets (~1M cells), `Sketch_parallel/` provides an alternative,
more memory-conservative route that geometrically sketches the data down to
50k representative cells, integrates the sketch, then projects the result
back onto the full dataset. Recommended chain, each stage reading the
previous stage's output:

1. `Sketch_parallel/main_pipeline.R` -- sketches `seurat_obj_pre_integration.rds`,
   hands off the sketch subset to `integrate_sketch.R` (run as a sub-process
   so the 1M-cell object isn't duplicated across forked workers), then
   projects and saves `seurat_obj_integrated_sketch.rds`.
   *(Alternative: `Sketch_parallel/FastMN.R` integrates via FastMNN instead of
   RPCA and saves `seurat_obj_mnn_integrated_1M.rds`; `cr_filtered/FastMN.R`
   is the same script for a CellRanger-filtered input.)*
2. `Sketch_parallel/Integrated_object_FastMNN/processed_reads/processing.R` --
   multi-resolution clustering + marker detection on the FastMNN-integrated
   object, saves `seurat_obj_clustered_1M.rds`.
3. `annotate.R` / `annotate_wu_suppiah.R` / `heatmaps.R` / `read_clustered.R` /
   `seurat_dimplot.R` -- annotation comparison plots, composition heatmaps,
   and per-resolution marker/heatmap generation from the clustered object.
4. `split_lab.R` + `lab_subsets/wu_lab.R` and the various `*_plot.R` /
   `diag.R` scripts -- per-lab subsetting, resolution-vs-annotation
   validation (ARI), and diagnostic/QC plots.

> **Before running:** `Sketch_parallel/*/FastMN.R` and
> `Sketch_parallel/cr_filtered/FastMN.R` call `ProjectIntegration()` with
> `reduction=`/`reduction.name=` arguments, while the RPCA/sketch scripts
> (`MPNST_RPCA_Integration.R`, `Sketch.R`, `main_pipeline.R`) use
> `integration.reduction=`/`new.reduction=` for the same call. Both forms
> have appeared in different Seurat versions -- verify which one matches
> your installed Seurat/SeuratWrappers version before running FastMN.R.

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

| File | Produced by | Description |
|------|-------------|-------------|
| `seurat_obj_pre_integration.rds` | MPNST.R | Post-QC Seurat object ready for integration |
| `01_QC_Pre_Scatter_Uniform.png` | MPNST.R | Pre-filter QC scatter plots |
| `01_QC_Pre_Filtering_Stats.csv` | MPNST.R | Pre-filter cell count statistics per batch/lab |
| `02_QC_Post_Scatter_Uniform.png` | MPNST.R | Post-filter QC scatter plots |
| `02_QC_Post_Filtering_Stats.csv` | MPNST.R | Post-filter cell count statistics per batch/lab |
| `03_UMAP_Unintegrated.png` | MPNST_RPCA_Integration.R | UMAP of unintegrated data |
| `04_UMAP_Integrated_RPCA.png` | MPNST_RPCA_Integration.R | UMAP after RPCA integration |
| `05_Integration_Comparison.png` | MPNST_RPCA_Integration.R | Side-by-side unintegrated vs. integrated UMAP |
| `seurat_obj_integrated_rpca.rds` | MPNST_RPCA_Integration.R | Final RPCA-integrated Seurat object |

> **Note:** the filenames above were corrected to match what the scripts
> actually write; a previous revision of this README listed
> `04_UMAP_Integrated.png`, `pre_filter_summary.csv`, and
> `post_filter_summary.csv`, none of which the code produces.

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

## Snakemake Workflow

Rather than running each R script or `.slurm` job by hand, you can drive the
pipeline through Snakemake, which tracks which stages are up to date and
resumes where it left off.

### 1. Install Snakemake (in addition to the R packages below)
```bash
pip install snakemake "pulp<2.8"
```

### 2. Configure the run
Edit `config/config.yaml`:
```yaml
input_rds: aggregated_cellbender_filtered.rds
integration: rpca   # rpca | sketch | sketch_parallel | fastmnn
run_clustering: false
```

### 3. Run
```bash
# Preview the jobs Snakemake will run
snakemake --cores 1 -n

# Run using R/Seurat already on PATH
snakemake --cores 4

# Or let Snakemake manage an isolated conda environment (envs/environment.yaml)
snakemake --cores 4 --use-conda
```

The `integration` option selects which batch-integration strategy feeds off
of `seurat_obj_pre_integration.rds` (see "Extended Pipeline" above for what
each one does). Setting `run_clustering: true` (only valid with
`integration: fastmnn`) additionally runs the multi-resolution clustering
and annotation-comparison stage under
`Sketch_parallel/Integrated_object_FastMNN/processed_reads/`.

> **Scale note:** this wraps the real pipeline, including the 30-40 hour
> RPCA integration step -- Snakemake changes how you launch and resume the
> work, not how long the underlying computation takes.

## Continuous Integration

`.github/workflows/ci.yml` runs on every push and pull request:

- **syntax-check** -- parses every `.R` file with base R and validates every
  `.slurm` script as bash, with no package installation required.
- **snakemake-dry-run** -- validates the Snakemake DAG resolves correctly
  for all four integration modes plus the clustering branch, and that
  invalid config combinations are rejected.
- **smoke-test** -- installs R + Seurat and runs the core QC + RPCA
  integration stages end to end against a tiny synthetic Seurat object
  (`tests/make_test_fixture.R`), then checks that every expected output
  file was produced. This does not (and cannot) exercise the full pipeline
  against real ~1M-cell data or the FastMNN/BPCells branches, which depend
  on packages and data not practical to install in CI; those are covered
  by the syntax and DAG-validation jobs instead.

## Author

**Vipin Menon**
Post doctoral Fellow, Computational Biology (Genome Editing), Wei Lei Lab
University of Maryland, Baltimore
University of Maryland-Institute of health computing 
GitHub: [vipinmenon1989](https://github.com/vipinmenon1989)
