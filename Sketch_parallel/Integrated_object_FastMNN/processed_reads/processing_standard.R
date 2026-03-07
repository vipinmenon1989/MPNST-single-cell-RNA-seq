#' @title Perform MNN-integrated Seurat analysis for clustering and marker detection.
#' @description This script loads a pre-integrated Seurat object, performs multi-resolution clustering,
#'   generates UMAP visualizations, and identifies marker genes for each cluster at various resolutions.
#' @details
#'   The pipeline includes:
#'   1. Setting memory limits for future globals to accommodate large datasets.
#'   2. Loading a previously saved Seurat object that has undergone MNN integration.
#'   3. Calculating k-nearest neighbors and performing multi-resolution Leiden/Louvain clustering
#'      based on the integrated reduction.
#'   4. Iterating through a range of clustering resolutions to:
#'      - Generate UMAP plots for each resolution, highlighting the identified clusters.
#'      - Identify differentially expressed marker genes for each cluster within that resolution.
#'   5. Saving a combined grid plot of all generated UMAPs and the final processed Seurat object.
#'
#'   Intermediate outputs include UMAP plots for each resolution (PNG), marker gene tables
#'   for each resolution (CSV), and a combined UMAP grid plot (PDF).
#' @param None This is a standalone script; all inputs are either file-based or configured internally.
#' @return
#'   This script does not return an R object directly but produces several output files:
#'   - `umap_outputs/UMAP_res_*.png`: Individual UMAP plots for each clustering resolution.
#'   - `umap_outputs/UMAP_Grid_All_Res.pdf`: A single PDF file containing a grid of all UMAP plots.
#'   - `marker_tables/markers_res_*.csv`: CSV files detailing marker genes for each cluster at each resolution.
#'   - `seurat_obj_clustered_1M.rds`: The final Seurat object, updated with clustering results, saved as an RDS file.
#' @seealso
#'   [Seurat::FindNeighbors()], [Seurat::FindClusters()], [Seurat::DimPlot()], [Seurat::FindAllMarkers()]
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import data.table
#' @keywords seurat, single-cell, RNA-seq, clustering, markers, integration, UMAP
#' @author Your Name
NULL

# 1. SET MEMORY LIMITS
options(future.globals.maxSize = 100 * 1024^3)

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)

# 2. Load Data
print("--- LOADING INTEGRATED OBJECT ---")
seurat_obj <- readRDS("seurat_obj_mnn_integrated_1M.rds")

# 3. Clustering
print("--- CALCULATING NEIGHBORS ---")
# Ensure the graph is named explicitly
seurat_obj <- FindNeighbors(
  seurat_obj, 
  reduction = "integrated.mnn.full", 
  dims = 1:30,
  graph.name = "integrated_snn" 
)

print("--- RUNNING MULTI-RESOLUTION CLUSTERING ---")
res_range <- seq(0, 1, by = 0.1)

seurat_obj <- FindClusters(
  seurat_obj, 
  graph.name = "integrated_snn",
  resolution = res_range
)

# 4. Visualization & Marker Detection
print("--- STARTING PLOT AND MARKER LOOP ---")
dir.create("umap_outputs", showWarnings = FALSE)
dir.create("marker_tables", showWarnings = FALSE)

plot_list <- list()

for(res in res_range) {
  
  # Construct column name
  res_col <- paste0("integrated_snn_res.", res)
  
  # --- CRITICAL FIX IS HERE ---
  # Use @meta.data[[ ]] to get a Vector, not a Data Frame
  if (!res_col %in% colnames(seurat_obj@meta.data)) {
    print(paste("SKIPPING:", res_col, "- not found in metadata"))
    next
  }
  
  # Assign Identity safely
  Idents(seurat_obj) <- seurat_obj@meta.data[[res_col]]
  
  print(paste0(">>> Processing Resolution: ", res))

  # A. Generate UMAP
  p <- DimPlot(seurat_obj, reduction = "umap.mnn", label = TRUE, raster = FALSE) + 
    ggtitle(paste("Resolution:", res)) +
    theme(legend.position = "none")
  
  ggsave(paste0("umap_outputs/UMAP_res_", res, ".png"), p, width = 10, height = 10)
  plot_list[[as.character(res)]] <- p
  
  # B. Find Markers
  print(paste0("   Finding markers for Res ", res, "..."))
  
  # Added min.cells.group to prevent crash on tiny singleton clusters
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    min.cells.group = 5,  # Skip tiny clusters (prevent errors)
    verbose = FALSE
  )
  
  write.csv(
    markers, 
    file = paste0("marker_tables/markers_res_", res, ".csv"), 
    row.names = FALSE
  )
  
  rm(markers)
  gc()
}

# 5. Save Grid Plot
print("--- SAVING GRID PLOT ---")
grid_plot <- wrap_plots(plot_list, ncol = 3)
ggsave("umap_outputs/UMAP_Grid_All_Res.pdf", grid_plot, width = 25, height = 35, limitsize = FALSE)



# 6. Save Final Object
saveRDS(seurat_obj, "seurat_obj_clustered_1M.rds")
print("--- PIPELINE SUCCESSFUL ---")
