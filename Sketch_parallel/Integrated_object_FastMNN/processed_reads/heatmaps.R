library(Seurat)
library(ggplot2)
library(reshape2)

# 1. INITIALIZE VARIABLES (Fixes 'object not found' errors)
setwd("/local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/Integrated_object_FastMNN/processed_reads/")
seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

resolutions <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
annot <- "wu_annot"

# 2. FILTER DATA
# Remove unknowns so they don't appear in the heatmap rows
cells_to_keep <- which(!seurat_obj[[annot, drop = TRUE]] %in% c("unknown", "Unknown", NA))
seurat_annot_filtered <- seurat_obj[, cells_to_keep]

# 3. GENERATE HEATMAPS
pdf(paste0("Cluster_Composition_Final_Heatmaps_", annot, ".pdf"), width = 16, height = 9)

for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  
  if (res_col %in% colnames(seurat_obj@meta.data)) {
    
    # Calculate contingency table
    ct <- as.data.frame(table(seurat_annot_filtered@meta.data[[res_col]], 
                              seurat_annot_filtered@meta.data[[annot]]))
    colnames(ct) <- c("Cluster", "Annotation", "Count")
    
    # Normalize by Cluster (Column)
    ct_wide <- dcast(ct, Annotation ~ Cluster, value.var = "Count", fill = 0)
    ct_mat <- as.matrix(ct_wide[,-1])
    rownames(ct_mat) <- ct_wide[,1]
    
    ct_prop <- sweep(ct_mat, 2, colSums(ct_mat), "/")
    ct_prop[is.nan(ct_prop)] <- 0 
    
    ct_plot <- melt(ct_prop)
    colnames(ct_plot) <- c("Annotation", "Cluster", "Proportion")
    
    # ENSURE EVERY CLUSTER IS LABELED: Treat Cluster as a factor in order
    ct_plot$Cluster <- factor(ct_plot$Cluster, levels = sort(as.numeric(unique(as.character(ct_plot$Cluster)))))

    # PLOT WITH BORDERS
    p_heat <- ggplot(ct_plot, aes(x = Cluster, y = Annotation, fill = Proportion)) +
              geom_tile(color = "black", linewidth = 0.2) + # Explicit black grid borders
              scale_fill_gradientn(colors = c("white", "yellow", "orange", "red", "darkred"), name = "Prop") +
              theme_bw() +
              labs(title = paste("Composition Analysis (Res:", res, ")"),
                   subtitle = "Grid borders added | Unknowns excluded | All clusters labeled",
                   x = "Cluster ID", y = "Cell Type") +
              theme(
                # Force every X-axis ID to be printed vertically
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
                axis.text.y = element_text(size = 10, face = "bold"),
                panel.grid = element_blank(),
                axis.ticks = element_line(color = "black")
              )

    print(p_heat)
  }
}
dev.off()