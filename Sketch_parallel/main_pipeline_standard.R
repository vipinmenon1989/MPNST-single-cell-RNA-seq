#' @file sketch_integration_pipeline.R
#' @title Seurat Sketch Integration Workflow for Large Datasets
#' @description This script orchestrates a Seurat integration pipeline specifically designed for very large datasets
#'   using the SketchData approach. It outlines the process of loading initial data, generating a sketch,
#'   performing integration on the slimmed-down sketch in an isolated process, and then projecting
#'   the integrated results back to the full dataset for UMAP visualization and downstream analysis.
#'
#' @details The pipeline proceeds through several critical stages:
#'   \enumerate{
#'     \item \strong{Setup}: Initializes the `future` plan for potential parallel processing.
#'     \item \strong{Load & Sketch}: Reads a pre-integration Seurat object from disk and creates a sketch
#'       (e.g., using LeverageScore) to efficiently represent the full dataset.
#'     \item \strong{Prepare Sketch Assay}: Standard Seurat preprocessing steps (FindVariableFeatures, ScaleData, RunPCA)
#'       are applied to the newly created sketch assay.
#'     \item \strong{Hand-off (Nuclear Option)}: To manage memory and computational resources, a highly slimmed-down
#'       Seurat object containing *only* the sketch assay and its PCA reduction is created and saved. This
#'       object is then passed to an external R process (`integrate_sketch.R`) for isolated integration.
#'       This step drastically reduces the memory footprint during the integration phase.
#'     \item \strong{Re-import and Merge}: After the external `integrate_sketch.R` script completes its
#'       task and saves the integrated sketch object, this script re-imports the integrated results.
#'       The integration reduction from the sketch is then seamlessly merged back into the primary (full)
#'       Seurat object. Temporary files are removed.
#'     \item \strong{Project & UMAP}: The integration performed on the sketch is projected onto the full
#'       RNA assay of the original Seurat object. Subsequently, UMAP dimensionality reduction is applied
#'       using the projected integrated data, providing a unified embedding for the entire dataset.
#'     \item \strong{Save Final}: The final, fully integrated Seurat object, now containing the projected
#'       integration and UMAP embeddings, is saved to disk for subsequent analysis.
#'   }
#'   This methodology is particularly beneficial for large-scale single-cell datasets where direct
#'   integration of the complete dataset might be computationally prohibitive or exceed available memory.
#'
#' @seealso \code{\link[Seurat]{SketchData}}, \code{\link[Seurat]{ProjectData}}, \code{\link[Seurat]{ProjectIntegration}}
#' @import Seurat
#' @import ggplot2
#' @import future
#' @import dplyr
#' @keywords Seurat Sketch Integration LargeData UMAP DimensionalityReduction
#'
library(Seurat)
library(ggplot2)
library(future)
library(dplyr)

# 1. SETUP
plan("sequential")
options(future.globals.maxSize = 450 * 1024^3) 

# 2. LOAD & SKETCH
print("--- LOADING DATA ---")
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")

DefaultAssay(seurat_obj) <- "RNA"
print("--- SKETCHING ---")
seurat_obj <- SketchData(
  object = seurat_obj, 
  ncells = 50000, 
  method = "LeverageScore",  
  sketched.assay = "sketch"
)

# 3. PREPARE SKETCH ASSAY
DefaultAssay(seurat_obj) <- "sketch"
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# 4. THE HAND-OFF (NUCLEAR OPTION)
print("--- EXPORTING SUBSET FOR ISOLATED INTEGRATION ---")
sketch_subset <- subset(seurat_obj, cells = Cells(seurat_obj[["sketch"]]))

# CRITICAL FIX: Strip the massive 'RNA' assay. Keep ONLY 'sketch' and 'pca'.
# This reduces the file size from ~250GB to ~2GB.
print("--- SLIMMING DOWN OBJECT (DietSeurat) ---")
sketch_subset <- DietSeurat(
  sketch_subset, 
  assays = "sketch", 
  dimreducs = "pca"
)

# Verify size before saving (Optional check for your logs)
print(format(object.size(sketch_subset), units = "GB"))

# Save with version 3
saveRDS(sketch_subset, "temp_sketch_only.rds", version = 3, compress = FALSE)
rm(sketch_subset)
gc()

print("--- LAUNCHING EXTERNAL R PROCESS (integrate_sketch.R) ---")
system("Rscript integrate_sketch.R")

# 5. RE-IMPORT AND MERGE
print("--- MERGING INTEGRATED RESULTS ---")
if(!file.exists("temp_sketch_integrated.rds")) {
  stop("Integration script failed to produce output.")
}
integrated_subset <- readRDS("temp_sketch_integrated.rds")
seurat_obj[["integrated.rpca"]] <- integrated_subset[["integrated.rpca"]]

# Cleanup temp files
file.remove("temp_sketch_only.rds")
file.remove("temp_sketch_integrated.rds")
rm(integrated_subset)
gc()

# 6. PROJECT & UMAP
print("--- PROJECTING TO FULL DATASET ---")
# Keep your current projection logic
seurat_obj <- ProjectData(
  object = seurat_obj,
  assay = "RNA",
  full.reduction = "pca",       
  sketched.assay = "sketch",
  sketched.reduction = "pca",    
  nmisc = 20
)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- ProjectIntegration(
  object = seurat_obj, 
  sketched.assay = "sketch", 
  integration.reduction = "integrated.rpca", 
  new.reduction = "integrated.rpca.full"
)

print("--- RUNNING UMAP ---")
seurat_obj <- RunUMAP(
  seurat_obj, 
  reduction = "integrated.rpca.full", 
  dims = 1:30, 
  reduction.name = "umap.full", 
  verbose = FALSE
)

# 7. SAVE FINAL
print("--- SAVING FINAL OBJECT ---")
saveRDS(seurat_obj, "seurat_obj_integrated_sketch.rds", version = 3, compress = FALSE)
print("--- PIPELINE COMPLETE ---")
