#' @title Compare Seurat Clustering Resolutions with WU Annotation
#'
#' @description
#' This script loads a pre-processed Seurat object, iterates through specified
#' clustering resolutions, and generates UMAP plots comparing each resolution's
#' clusters with a "wu_annot" reference annotation.
#'
#' @details
#' The script expects a Seurat object file named 'seurat_obj_clustered_1M.rds'
#' to be present in the working directory or provided via a full path.
#' This Seurat object should contain:
#' \itemize{
#'   \item A 'umap.mnn' dimensional reduction.
#'   \item Metadata columns for various 'integrated_snn_res.*' resolutions
#'         (e.g., 'integrated_snn_res.0.1', 'integrated_snn_res.0.2', etc.).
#'   \item A 'wu_annot' metadata column for reference annotation.
#' }
#' For the annotation plot, cells with "unknown" or NA values in 'wu_annot'
#' are explicitly excluded to create a cleaner reference. All UMAP plots are
#' rasterized to efficiently handle large cell numbers (e.g., 1M cells)
#' within the PDF output, preventing potential rendering issues with vector
#' graphics. Memory cleanup (`gc()`) is performed within the loop to manage
#' resources for large datasets.
#'
#' @param seurat_obj_clustered_1M.rds A character string specifying the path
#'   to the input Seurat object file. This is an implicit input to the script.
#'   The default expectation is a file named "seurat_obj_clustered_1M.rds"
#'   in the current working directory.
#'
#' @return
#' `resolutions_comparison_final.pdf`: A PDF file containing a series of
#' paired UMAP plots. Each pair displays:
#' \enumerate{
#'   \item The full clustering for a given resolution (including all cells).
#'   \item The filtered WU annotation plot (excluding "unknown" and NA values).
#' }
#' The PDF is saved to the current working directory.
#'
#' @examples
#' # This script is designed to be run directly.
#' # Before running, ensure 'seurat_obj_clustered_1M.rds' is available
#' # in the working directory or adjust the `readRDS()` path accordingly.
#' #
#' # Example usage (assuming the script is saved as 'compare_resolutions.R'):
#' # source("compare_resolutions.R")
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
library(Seurat)
library(ggplot2)
library(patchwork)

# 1. Load the object
seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

# 2. Define resolutions
resolutions <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")

# 3. Create the filtered subset for the annotation plot
# This logic removes:
# - Any cell where wu_annot is literally the string "unknown"
# - Any cell where wu_annot is NA (missing data)
seurat_annot_only <- subset(seurat_obj, 
                            subset = wu_annot != "unknown" & !is.na(wu_annot))

# 4. Open PDF with rasterization for 1M cells
pdf("resolutions_comparison_final.pdf", width = 16, height = 8)

for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  
  if (res_col %in% colnames(seurat_obj@meta.data)) {
    
    # Plot 1: Full Clustering (All cells, including NA/Unknown)
    # We use raster=TRUE because 1M vector points will break your PDF viewer
    p1 <- DimPlot(seurat_obj, 
                  reduction = "umap.mnn", 
                  group.by = res_col, 
                  label = TRUE, 
                  label.size = 3,
                  raster = TRUE) + 
      ggtitle(paste("Full Resolution:", res)) +
      theme(legend.position = "none")

    # Plot 2: Filtered Annotation (No Unknown, No NA)
    p2 <- DimPlot(seurat_annot_only, 
                  reduction = "umap.mnn", 
                  group.by = "wu_annot", 
                  label = TRUE, 
                  label.size = 3,
                  raster = TRUE) + 
      ggtitle("WU Annotation (Filtered: No Unknown/NA)")

    # Print to PDF
    print(p1 + p2)
    
    # Memory cleanup inside the loop
    gc()
    message(paste("Processed resolution:", res))
  }
}

dev.off()
