library(Seurat)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("/local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/Integrated_object_FastMNN/processed_reads/")

# Load object
seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

# Configuration
resolutions <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
annot <- "wu_annot"

# Filter for Annotation UMAP (removing Unknown/NA)
cells_to_keep <- which(!seurat_obj[[annot, drop = TRUE]] %in% c("unknown", "Unknown", NA))
seurat_annot_filtered <- seurat_obj[, cells_to_keep]

# Generate UMAP-only PDF
pdf(paste0("UMAP_Comparison_Clean_", annot, ".pdf"), width = 18, height = 8)

for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  
  if (res_col %in% colnames(seurat_obj@meta.data)) {
    
    # Left Plot: Clusters (All cells)
    p1 <- DimPlot(seurat_obj, reduction = "umap.mnn", group.by = res_col, 
                  label = TRUE, raster = TRUE) + 
          ggtitle(paste("Resolution:", res)) + 
          theme(legend.position = "none")

    # Right Plot: Wu Annotations (Filtered)
    p2 <- DimPlot(seurat_annot_filtered, reduction = "umap.mnn", group.by = annot, 
                  label = TRUE, raster = TRUE) + 
          ggtitle(paste("Annotation:", annot))

    # Print only the two UMAPs side-by-side
    print(p1 + p2)
  }
}
dev.off()