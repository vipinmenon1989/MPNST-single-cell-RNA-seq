#' @title Seurat FastMNN Integration Pipeline with Geometric Sketching for Large Datasets
#' @description This script provides a complete pipeline for integrating large-scale
#'   single-cell RNA-seq datasets (e.g., 1M cells) using Seurat's FastMNN integration
#'   method, enhanced by geometric sketching to reduce computational burden.
#'   It is optimized for High-Performance Computing (HPC) environments by
#'   strategically managing memory and parallel processing.
#'
#'   The core strategy involves:
#'   \enumerate{
#'     \item Loading a pre-integration Seurat object.
#'     \item Creating a small (e.g., 50k cells) geometric sketch to represent the full dataset.
#'     \item Evicting the full object from RAM and performing intensive processing (variable features, scaling, PCA)
#'           and FastMNN integration on the sketch subset.
#'     \item Reloading the full object and projecting the integration results from the sketch
#'           onto the entire dataset, effectively integrating 1M cells based on the
#'           manifold learned from the 50k sketch.
#'     \item Performing final dimension reduction (UMAP) and saving the integrated object.
#'   }
#'   This approach drastically reduces memory requirements and computational time
#'   for large datasets compared to direct integration of all cells.
#'
#' @details This script uses `future` for parallelization, `Seurat` for single-cell
#'   data analysis, `SeuratWrappers` for `FastMNNIntegration`, and `batchelor`
#'   as the underlying engine for FastMNN.
#'
#' @section HPC Considerations:
#'   The script initiates with a sequential `future` plan for safe loading of large
#'   objects, then switches to `multicore` for parallel FastMNN integration once
#'   the large object is offloaded from RAM. Memory limits for global variables
#'   are set via `options(future.globals.maxSize)`.
#'
#' @section Input:
#'   \code{seurat_obj_pre_integration.rds}: An RDS file containing a Seurat object
#'   with multiple assays/layers representing different batches or samples,
#'   prepared for integration.
#'
#' @section Output:
#'   \code{seurat_obj_mnn_integrated_1M.rds}: An RDS file containing the final
#'   integrated Seurat object with a new "integrated.mnn.full" reduction and a
#'   "umap.mnn" UMAP reduction.
#'
#' @author [Your Name Here]
#' @seealso \code{\link[Seurat]{SketchData}}, \code{\link[SeuratWrappers]{FastMNNIntegration}},
#'   \code{\link[Seurat]{IntegrateLayers}}, \code{\link[Seurat]{ProjectIntegration}}
#' @concept Seurat
#' @concept FastMNN
#' @concept single-cell RNA-seq
#' @concept large data integration
#' @concept geometric sketching
#' @concept HPC
#'
#' @import Seurat
#' @import SeuratWrappers
#' @import batchelor
#' @import ggplot2
#' @import patchwork
#' @import future
#' @import dplyr
NULL

library(Seurat)
library(SeuratWrappers) # Required for FastMNNIntegration
library(batchelor)      # The underlying engine for FastMNN
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

# 1. HPC SETUP
# -------------------------------------------------------------------------
# Using 'multicore' for Linux HPC; use 'multisession' if on Windows.
plan("sequential")
options(future.globals.maxSize = 400 * 1024^3) # 400GB limit for globals

print("--- HPC MODE: INITIALIZED ---")

# 2. LOAD DATA
# -------------------------------------------------------------------------
file_path <- "seurat_obj_pre_integration.rds"

if (!file.exists(file_path)) {
  stop("Input file not found.")
}

seurat_obj <- readRDS(file_path)
print(paste("--- DATA LOADED: ", ncol(seurat_obj), " cells ---"))

# 3. GEOMETRIC SKETCHING
# -------------------------------------------------------------------------
# We sketch 50k cells to represent the 1M population.
# This is the "Blunt Truth" for speed: 1M cell integration is unnecessary
# if a 50k sketch captures the manifold.
print("--- RUNNING GEOMETRIC SKETCH (50k cells) ---")
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- SketchData(
  object = seurat_obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# 4. MEMORY EVICTION & SUBSETTING
# -------------------------------------------------------------------------
print("--- SUBSETTING SKETCH & EVICTING FULL OBJECT ---")
sketch_subset <- subset(seurat_obj, cells = Cells(seurat_obj[["sketch"]]))
DefaultAssay(sketch_subset) <- "sketch"

# Save the full object to disk and WIPE from RAM
saveRDS(seurat_obj, "temp_seurat_full.rds", compress = FALSE)
rm(seurat_obj)
gc() # Deep garbage collection

# 5. PROCESS SKETCH
# -------------------------------------------------------------------------
print("--- PROCESSING SKETCH SUBSET ---")
sketch_subset <- FindVariableFeatures(sketch_subset, verbose = FALSE)
sketch_subset <- ScaleData(sketch_subset, verbose = FALSE)
sketch_subset <- RunPCA(sketch_subset, verbose = FALSE)

# 6. PARALLEL FastMNN INTEGRATION
# -------------------------------------------------------------------------
# Switching to multicore NOW because the environment is empty of the 1M cell object.
print("--- HPC MODE: MULTICORE (8 WORKERS) ---")
plan("multicore", workers = 8)

print("--- INTEGRATING LAYERS VIA FastMNN ---")
# FastMNN integrates in the PCA space; it is faster and more stable than RPCA for high-K batches.
sketch_subset <- IntegrateLayers(
  object = sketch_subset,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = TRUE
)

# 7. RE-LINKING & PROJECTION
# -------------------------------------------------------------------------
plan("sequential")
print("--- RELOADING MAIN OBJECT & PROJECTING ---")
seurat_obj <- readRDS("temp_seurat_full.rds")

# Transfer integration from sketch to full
seurat_obj[["integrated.mnn"]] <- sketch_subset[["integrated.mnn"]]
rm(sketch_subset)
gc()

seurat_obj <- ProjectIntegration(
  object = seurat_obj,
  sketched.assay = "sketch",
  reduction = "integrated.mnn",           # Changed from integration.reduction
  reduction.name = "integrated.mnn.full"  # Changed from new.reduction
)

# 8. FINAL DIM REDUCTION & SAVE
# -------------------------------------------------------------------------ls
print("--- FINAL UMAP FOR 1M CELLS ---")
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "integrated.mnn.full",
  dims = 1:30,
  reduction.name = "umap.mnn",
  verbose = FALSE
)

print("--- SAVING FINAL INTEGRATED OBJECT ---")
saveRDS(seurat_obj, "seurat_obj_mnn_integrated_1M.rds", compress = FALSE)
file.remove("temp_seurat_full.rds")

print("--- PIPELINE COMPLETED ---")
