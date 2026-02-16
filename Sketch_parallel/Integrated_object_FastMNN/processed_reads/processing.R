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