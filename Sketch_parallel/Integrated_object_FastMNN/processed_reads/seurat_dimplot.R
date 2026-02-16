# 1. SET MEMORY LIMITS
options(future.globals.maxSize = 500 * 1024^3)

library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

# 2. LOAD DATA
print("--- LOADING OBJECT ---")
 seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

# CRITICAL: Do NOT run JoinLayers(seurat_obj) here. It will crash.
DefaultAssay(seurat_obj) <- "RNA"

# Create output directories
dir.create("marker_tables_final", showWarnings = FALSE)
dir.create("heatmap_outputs_final", showWarnings = FALSE)

# 3. DEFINE RESOLUTIONS
res_range <- seq(0.1, 1, by = 0.1)

print("--- STARTING SUBSET-FIRST PIPELINE ---")

for (res in res_range) {
  
  res_col <- paste0("integrated_snn_res.", res)
  print(paste(">>> PROCESSING RESOLUTION:", res))
  
  # A. VALIDATE METADATA
  if (!res_col %in% colnames(seurat_obj@meta.data)) {
    print(paste("   Skipping:", res_col, "- not found."))
    next
  }
  
  # Set Identity to current resolution
  Idents(seurat_obj) <- seurat_obj@meta.data[[res_col]]
  
  # B. SMART SUBSET (THE FIX)
  # We extract max 500 cells per cluster.
  # This reduces the matrix size from 1,000,000 to ~15,000 (if 30 clusters).
  # This is statistically sufficient for marker detection and safe for memory.
  print("   Subsetting 500 cells/cluster to bypass integer limit...")
  
  cells_to_keep <- unlist(lapply(levels(Idents(seurat_obj)), function(x) {
    cells_in_cluster <- WhichCells(seurat_obj, idents = x)
    if (length(cells_in_cluster) > 500) {
      sample(cells_in_cluster, 500)
    } else {
      cells_in_cluster
    }
  }))
  
  # Create a lightweight temporary object
  sub_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  # C. JOIN LAYERS ON SUBSET
  # This works because sub_obj is small.
  # It unifies the split layers into one 'counts' matrix for FindAllMarkers.
  sub_obj <- JoinLayers(sub_obj)
  
  # D. FIND MARKERS
  print("   Calculating markers on subset...")
  
  tryCatch({
    # No need for max.cells.per.ident here because sub_obj is already downsampled
    markers <- FindAllMarkers(
      sub_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = FALSE
    )
    
    if (nrow(markers) == 0) {
      print("   WARNING: No markers found. Skipping.")
      next
    }
    
    # Save Markers
    write.csv(
      markers, 
      file = paste0("marker_tables_final/markers_res_", res, ".csv"),
      row.names = FALSE
    )
    
    # E. PREPARE HEATMAP DATA
    if (!"cluster" %in% colnames(markers)) next
    
    top5_genes <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 5, order_by = avg_log2FC) %>%
      pull(gene) %>%
      unique()
      
    if (length(top5_genes) == 0) next
    
    # F. GENERATE HEATMAP
    # We downsample again to 50 cells/cluster for the visual plot (cleaner)
    print("   Generating Heatmap...")
    heatmap_cells <- unlist(lapply(levels(Idents(sub_obj)), function(x) {
      cells <- WhichCells(sub_obj, idents = x)
      if (length(cells) > 50) sample(cells, 50) else cells
    }))
    
    heatmap_obj <- subset(sub_obj, cells = heatmap_cells)
    heatmap_obj <- ScaleData(heatmap_obj, features = top5_genes, verbose = FALSE)
    
    p <- DoHeatmap(
      heatmap_obj, 
      features = top5_genes, 
      size = 3, 
      angle = 90
    ) + 
    ggtitle(paste("Resolution:", res)) +
    theme(axis.text.y = element_text(size = 5))

    calc_height <- max(8, length(top5_genes) * 0.15)
    
    ggsave(
      filename = paste0("heatmap_outputs_final/Heatmap_res_", res, ".pdf"),
      plot = p, 
      width = 15, 
      height = calc_height,
      limitsize = FALSE
    )
    
    print(paste("   SUCCESS: Saved Res", res))
    
  }, error = function(e) {
    print(paste("   ERROR at Res", res, ":", e$message))
  })
  
  # Cleanup to free RAM
  rm(sub_obj, heatmap_obj, markers)
  gc()
}

print("--- PIPELINE FINISHED ---")