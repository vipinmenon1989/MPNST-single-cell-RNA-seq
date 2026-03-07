#' @file seurat_sketch_integration_workflow.R
#' @title Seurat Integration Workflow for Large Datasets using Geometric Sketching
#' @description This script orchestrates a Seurat object integration workflow tailored
#'   for large datasets using geometric sketching. It leverages the `future` package
#'   for efficient parallel processing on HPC environments, particularly during
#'   the integration phase.
#' @details The workflow encompasses several key steps:
#'   1.  Initial setup for HPC mode using `future`.
#'   2.  Loading of a pre-integration Seurat object.
#'   3.  Generating a geometric sketch of the dataset for computational efficiency.
#'   4.  Optimized, parallel integration of the sketched data using RPCA,
#'       detaching the sketch from the main object to manage memory.
#'   5.  Projecting the learned integration to the full, original dataset.
#'   6.  Running UMAP on the fully integrated data and generating a visualization.
#'   This approach is designed to handle Seurat objects that are too large
#'   to be processed entirely in memory or sequentially.
#' @author [Your Name/Organization Here]
#' @date [Current Date, e.g., 2023-10-27]
#' @references Seurat V5 for large dataset integration and sketching methods.
#'
library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

#' @section 1. HPC Setup (Initial):
#' Initializes the `future` package for parallel processing.
#' Sets the plan to "sequential" initially to safely load large data
#' and configures `future` options to handle large global sizes.
#' This allows for safe data loading before enabling parallelization.
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3)

print("--- HPC MODE: INITIALIZED ---")

#' @section 2. Load Data:
#' Loads the pre-integration Seurat object from an RDS file.
#'
#' @var file_path Character string: Path to the input Seurat object RDS file.
file_path <- "seurat_obj_pre_integration.rds"
#' @var batch_column Character string: Specifies the metadata column in the Seurat object
#'   that identifies different batches or samples for integration and visualization.
batch_column <- "orig.ident"

if (!file.exists(file_path)) {
  stop("Input file not found. Ensure Part 1 ran successfully.")
}

#' Reads the Seurat object from the specified `file_path`.
#' The loaded object is stored in the `seurat_obj` variable.
seurat_obj <- readRDS(file_path)
print("--- DATA LOADED ---")

#' @section 3. Sketch the Dataset:
#' Performs geometric sketching on the Seurat object to create a smaller,
#' representative subset of cells. This sketch is then used for efficient
#' downstream integration.
print("--- RUNNING GEOMETRIC SKETCH ---")
#' Sets the default assay of the Seurat object to "RNA" prior to sketching.
DefaultAssay(seurat_obj) <- "RNA"

#' Applies the `SketchData` function to the Seurat object.
#' Creates a new assay named "sketch" containing 50,000 cells selected
#' using the "LeverageScore" method. This reduces the computational burden
#' for subsequent integration steps.
#' @note Warnings related to 'nsketch' can typically be ignored.
seurat_obj <- SketchData(
  object = seurat_obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

print("--- SKETCH COMPLETE ---")


#' @section 4. Optimized Integration (Detach & Parallelize):
#' This section performs integration of the sketched data. To manage memory
#' and leverage parallel processing, the sketch is temporarily detached
#' into a separate, lightweight Seurat object.
print("--- OPTIMIZATION: DETACHING SKETCH FOR PARALLEL PROCESSING ---")

#' @subsection A. Create a lightweight subset:
#' Creates `sketch_subset`, a new Seurat object containing only the cells
#' that were included in the "sketch" assay. This significantly reduces
#' the memory footprint when data is exported to parallel workers.
sketch_subset <- subset(seurat_obj, cells = Cells(seurat_obj[["sketch"]]))
#' Sets the default assay of the `sketch_subset` to "sketch".
DefaultAssay(sketch_subset) <- "sketch"

#' @subsection B. Pre-process the subset:
#' Performs standard Seurat preprocessing steps on the `sketch_subset`
#' in preparation for integration.
print("--- PROCESSING SKETCH SUBSET ---")
#' Finds variable features within the `sketch_subset`.
sketch_subset <- FindVariableFeatures(sketch_subset, verbose = FALSE)
#' Scales the data in the `sketch_subset`.
sketch_subset <- ScaleData(sketch_subset, verbose = FALSE)
#' Runs Principal Component Analysis (PCA) on the scaled `sketch_subset` data.
sketch_subset <- RunPCA(sketch_subset, verbose = FALSE)

#' @subsection C. Enable Parallelization:
#' Switches the `future` plan to "multicore" mode, enabling the use of
#' multiple CPU cores for the `IntegrateLayers` function.
print("--- HPC MODE: SWITCHING TO MULTICORE (8 WORKERS) ---")
plan("multicore", workers = 8)

print("--- INTEGRATING SKETCH LAYERS (RPCA - PARALLEL) ---")
#' Integrates the layers (batches) within the `sketch_subset` using RPCAIntegration.
#' This critical step now runs in parallel across 8 worker processes,
#' significantly reducing computation time. The new reduction is named "integrated.rpca".
sketch_subset <- IntegrateLayers(
  object = sketch_subset,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

#' @subsection D. Revert to Sequential & Reattach:
#' Reverts the `future` plan to sequential mode and transfers the
#' integration results from the `sketch_subset` back to the main `seurat_obj`.
plan("sequential")
print("--- HPC MODE: REVERTED TO SEQUENTIAL ---")

print("--- TRANSFERRING INTEGRATION TO MAIN OBJECT ---")
#' Copies the "integrated.rpca" reduction, computed on the sketch,
#' back into the `seurat_obj`.
seurat_obj[["integrated.rpca"]] <- sketch_subset[["integrated.rpca"]]

#' Cleans up the temporary `sketch_subset` object to free up system memory.
rm(sketch_subset)
gc()

print("--- SKETCH INTEGRATION COMPLETE ---")


#' @section 5. Project Integration:
#' Projects the integration, which was performed on the subset of sketched cells,
#' to the entire dataset. This maps all cells in the original object onto
#' the integrated low-dimensional space.
print("--- PROJECTING INTEGRATION TO FULL DATASET ---")
#' Sets the default assay of the main `seurat_obj` to "RNA" for projection.
DefaultAssay(seurat_obj) <- "RNA"

#' Applies the `ProjectIntegration` function.
#' This step uses the "sketch" assay and the "integrated.rpca" reduction
#' to create a new reduction, "integrated.rpca.full", that covers all cells
#' in the original `seurat_obj`.
seurat_obj <- ProjectIntegration(
  object = seurat_obj,
  sketched.assay = "sketch",
  integration.reduction = "integrated.rpca",
  new.reduction = "integrated.rpca.full"
)

print("--- PROJECTION COMPLETE ---")


#' @section 6. Run UMAP on Full Dataset:
#' Computes a Uniform Manifold Approximation and Projection (UMAP)
#' on the fully integrated and projected dataset for visualization.
print("--- RUNNING UMAP ON PROJECTED DATA ---")
#' Runs UMAP on the `seurat_obj` using the "integrated.rpca.full" reduction.
#' The result is stored as a new reduction named "umap.full".
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "integrated.rpca.full",
  dims = 1:30,
  reduction.name = "umap.full",
  verbose = FALSE
)

print("--- GENERATING INTEGRATED PLOT ---")
#' Generates a `DimPlot` of the UMAP results.
#' The plot displays the "umap.full" reduction, with cells grouped by
#' the `batch_column` to visualize integration quality.
p2 <- DimPlot(seurat_obj,
              reduction = "umap.full",
              group.by = batch_column,
              label = TRUE,
              label.size = 2,
              repel = TRUE,
              raster = FALSE) +
      NoLegend
