#' @file run_seurat_integration_pipeline.R
#' @title Seurat Integration Pipeline using Sketching and Projection (Sequential)
#' @description This script performs Seurat data integration for large datasets using a sketching and projection approach.
#'   It is specifically configured to run in sequential mode on a single core to avoid memory export errors
#'   often encountered with parallel processing on HPC systems for very large objects.
#'   The pipeline involves data sketching, integration of the sketch, projection of integration to the full dataset,
#'   UMAP dimensionality reduction, and saving the final integrated object.
#' @author N/A
#' @dependencies Seurat, ggplot2, patchwork, future, dplyr
#' @details This script is optimized for memory efficiency by forcing sequential processing via `plan("sequential")`.
#'   It uses the `SketchData` and `ProjectIntegration` functions from Seurat to handle large cell numbers efficiently.
#'   The output includes a UMAP plot of the integrated data and a saved Seurat object.
#' @note Ensure that the input file "seurat_obj_pre_integration.rds" exists in the current working directory.
#'   The `future.globals.maxSize` option is set high as a precaution, though sequential mode generally avoids these issues.
#'   Large objects are saved with `compress = FALSE` for faster I/O on HPC.

library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

# 1. HPC SETUP (SEQUENTIAL MODE)
#' @section 1. HPC Setup (Sequential Mode):
#'   Configures the `future` package for sequential processing to manage memory for large Seurat objects.
#'   This step is crucial for preventing common memory-related export errors on HPC.
#'   The global maximum size for future globals is also increased.
# -------------------------------------------------------------------------
# CRITICAL: We force sequential processing. 
# This bypasses the "304.77 GiB" export error completely.
plan("sequential") 

# We set the limit high just in case, but sequential mode shouldn't trigger it.
options(future.globals.maxSize = 500 * 1024^3)

print("--- HPC MODE: SEQUENTIAL (1 CORE) ---")

# 2. LOAD DATA
#' @section 2. Load Data:
#'   Loads the pre-integration Seurat object from an RDS file.
#'   A check is performed to ensure the input file exists before proceeding.
# -------------------------------------------------------------------------
file_path <- "seurat_obj_pre_integration.rds" 
batch_column <- "orig.ident"

if (!file.exists(file_path)) {
  stop("Input file not found. Ensure Part 1 ran successfully.")
}

seurat_obj <- readRDS(file_path)
print("--- DATA LOADED ---")

# 3. SKETCH THE DATASET
#' @section 3. Sketch the Dataset:
#'   Applies geometric sketching to the Seurat object using the Leverage Score method.
#'   This reduces the dataset size to a manageable number of cells (e.g., 50,000) for faster initial processing.
# -------------------------------------------------------------------------
print("--- RUNNING GEOMETRIC SKETCH (SEQUENTIAL) ---")
# 1 Core is sufficient here.

DefaultAssay(seurat_obj) <- "RNA"

# Note: Ignore 'nsketch' warnings.
seurat_obj <- SketchData(
  object = seurat_obj, 
  ncells = 50000, 
  method = "LeverageScore",  
  sketched.assay = "sketch"
)

print("--- SKETCH COMPLETE ---")


# 4. INTEGRATE THE SKETCH
#' @section 4. Integrate the Sketch:
#'   Performs initial processing (variable feature finding, scaling, PCA) and then integrates the sketched dataset.
#'   Integration is done using RPCAIntegration, which is suitable for batch correction.
# -------------------------------------------------------------------------
print("--- PROCESSING SKETCH LAYERS ---")
# We switch to the sketch assay (50k cells)
DefaultAssay(seurat_obj) <- "sketch"

# Basic processing on the sketch
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

print("--- INTEGRATING SKETCH LAYERS (RPCA - SEQUENTIAL) ---")
# On 1 core, integrating 50k cells takes about 30-60 mins.
seurat_obj <- IntegrateLayers(
  object = seurat_obj, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

print("--- SKETCH INTEGRATION COMPLETE ---")


# 5. PROJECT INTEGRATION
#' @section 5. Project Integration:
#'   Projects the integration results from the sketched dataset back onto the full dataset.
#'   This step efficiently maps the learned batch correction to all cells without re-running full integration.
# -------------------------------------------------------------------------
print("--- PROJECTING INTEGRATION TO FULL DATASET ---")
# This maps the 950k other cells onto the 50k anchors.
# This uses Matrix Algebra, which R handles efficiently on 1 core.

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- ProjectIntegration(
  object = seurat_obj, 
  sketched.assay = "sketch", 
  integration.reduction = "integrated.rpca", 
  new.reduction = "integrated.rpca.full"
)

print("--- PROJECTION COMPLETE ---")


# 6. RUN UMAP ON FULL DATASET
#' @section 6. Run UMAP on Full Dataset:
#'   Performs Uniform Manifold Approximation and Projection (UMAP) on the full, projected dataset.
#'   Generates and saves a UMAP plot grouped by the original batch column.
# -------------------------------------------------------------------------
print("--- RUNNING UMAP ON PROJECTED DATA ---")
# We use the full projected reduction 'integrated.rpca.full'
seurat_obj <- RunUMAP(
  seurat_obj, 
  reduction = "integrated.rpca.full", 
  dims = 1:30, 
  reduction.name = "umap.full", 
  verbose = FALSE
)

print("--- GENERATING INTEGRATED PLOT ---")
p2 <- DimPlot(seurat_obj, 
              reduction = "umap.full", 
              group.by = batch_column, 
              label = TRUE, 
              label.size = 2, 
              repel = TRUE, 
              raster = FALSE) + 
      NoLegend() + 
      ggtitle("Integrated (Sketch-Projected)") +
      theme(plot.title = element_text(hjust = 0.5))

ggsave("04_UMAP_Integrated_Sketch.png", plot = p2, width = 16, height = 12, dpi = 300)
print("--- INTEGRATED UMAP SAVED ---")


# 7. SAVE FINAL CHECKPOINT
#' @section 7. Save Final Checkpoint:
#'   Saves the final integrated Seurat object to an RDS file without compression.
#'   This creates a checkpoint for subsequent analysis steps.
# -------------------------------------------------------------------------
print("--- SAVING INTEGRATED OBJECT (NO COMPRESSION) ---")
# Always compress=FALSE for large objects on HPC
saveRDS(seurat_obj, "seurat_obj_integrated_sketch.rds", compress = FALSE)

print("--- PIPELINE COMPLETED ---")
