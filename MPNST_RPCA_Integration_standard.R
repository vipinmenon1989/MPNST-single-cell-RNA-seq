#' @file integrate_seurat_rpca.R
#' @title Seurat Object Integration with RPCA (Sequential Mode)
#' @description This script performs Seurat object integration using RPCA, optimized for large datasets on HPC environments by enforcing sequential processing.
#'   It loads a pre-integration Seurat object, runs unintegrated UMAP, performs RPCA integration,
#'   generates an integrated UMAP, and saves the final integrated object.
#' @details This script is designed to handle very large Seurat objects by limiting parallel processing
#'   to a single core during the `IntegrateLayers` step, preventing excessive RAM spikes.
#'   Estimated runtime is 30-40 hours for typical large datasets.
#' @author [Your Name/Organization Here]
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import future
library(Seurat)
library(ggplot2)
library(patchwork)
library(future)

#' 1. HPC SETUP (SEQUENTIAL MODE)
#' -------------------------------------------------------------------------
#' CRITICAL CHANGE: We use "sequential" to force 1 core.
#' This ensures RAM usage stays ~125GB and never spikes to 500GB.
plan("sequential")
options(future.globals.maxSize = 300 * 1024^3)
print("--- HPC MODE: SEQUENTIAL (1 CORE) ENABLED ---")
print("--- ESTIMATED RUNTIME: 30-40 HOURS ---")

#' 2. LOAD DATA
#' -------------------------------------------------------------------------
file_path <- "seurat_obj_pre_integration.rds"
batch_column <- "orig.ident"

if (!file.exists(file_path)) {
  stop("Input file not found. Ensure Part 1 ran successfully.")
}

seurat_obj <- readRDS(file_path)
print("--- DATA LOADED ---")

#' 3. UNINTEGRATED UMAP
#' -------------------------------------------------------------------------
print("--- RUNNING UMAP ON RAW PCA (UNINTEGRATED) ---")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated", verbose = FALSE)

print("--- GENERATING UNINTEGRATED PLOT ---")
p1 <- DimPlot(seurat_obj,
              reduction = "umap.unintegrated",
              group.by = batch_column,
              label = TRUE,
              label.size = 2,
              repel = TRUE,
              raster = FALSE) +
      NoLegend() +
      ggtitle("Unintegrated: Raw PCA") +
      theme(plot.title = element_text(hjust = 0.5))

ggsave("03_UMAP_Unintegrated.png", plot = p1, width = 16, height = 12, dpi = 300)
print("--- UNINTEGRATED UMAP SAVED ---")


#' 4. INTEGRATION (SEQUENTIAL RPCA)
#' -------------------------------------------------------------------------
print("--- RUNNING RPCA INTEGRATION (SEQUENTIAL) ---")
#' This is the heavy step. It will proceed 1 batch at a time.
#' Do not cancel if it looks "stuck" for an hour; it is just calculating.

seurat_obj <- IntegrateLayers(
  object = seurat_obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  k.weight = 100,
  verbose = TRUE # Keep verbose=TRUE to track progress in the log
)

print("--- INTEGRATION COMPLETE ---")


#' 5. INTEGRATED UMAP
#' -------------------------------------------------------------------------
print("--- RUNNING UMAP ON INTEGRATED DATA ---")
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", verbose = FALSE)

print("--- GENERATING INTEGRATED PLOT ---")
p2 <- DimPlot(seurat_obj,
              reduction = "umap.rpca",
              group.by = batch_column,
              label = TRUE,
              label.size = 2,
              repel = TRUE,
              raster = FALSE) +
      NoLegend() +
      ggtitle("Integrated (RPCA): Corrected") +
      theme(plot.title = element_text(hjust = 0.5))

ggsave("04_UMAP_Integrated_RPCA.png", plot = p2, width = 16, height = 12, dpi = 300)
print("--- INTEGRATED UMAP SAVED ---")


#' 6. SIDE-BY-SIDE COMPARISON
#' -------------------------------------------------------------------------
ggsave("05_Integration_Comparison.png", plot = p1 + p2, width = 24, height = 12, dpi = 300)


#' 7. SAVE FINAL CHECKPOINT
#' -------------------------------------------------------------------------
print("--- SAVING INTEGRATED OBJECT (NO COMPRESSION) ---")
#' Always compress=FALSE for large objects on HPC
saveRDS(seurat_obj, "seurat_obj_integrated_rpca.rds", compress = FALSE)

print("--- PIPELINE COMPLETED ---")
