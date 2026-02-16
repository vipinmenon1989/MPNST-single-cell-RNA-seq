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
