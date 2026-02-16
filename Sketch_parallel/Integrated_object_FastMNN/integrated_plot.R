library(Seurat)
library(ggplot2)

# 1. PREPARE LAB LIST
print("--- FETCHING LAB LIST ---")
temp_obj <- readRDS("seurat_obj_mnn_integrated_1M.rds")
all_labs <- unique(as.character(temp_obj$lab))
all_labs <- all_labs[!is.na(all_labs)]
rm(temp_obj); gc()

# 2. INTEGRATED LOOP
for (current_lab in all_labs) {
  print(paste0("--- INTEGRATED (CLEAN): ", current_lab, " ---"))
  
  # Load, Subset, and EVICT parent immediately
  xbj <- readRDS("seurat_obj_mnn_integrated_1M.rds")
  xbj_lab <- subset(xbj, subset = lab == current_lab)
  rm(xbj); gc()

  # Plotting (Using existing umap.mnn)
  for (annot in c("wu_annot", "suppiah_annot")) {
    if (annot %in% colnames(xbj_lab@meta.data)) {
      
      # --- THE FIX: REMOVE UNKNOWNS ---
      # Identify cells that are NOT 'unknown' or 'NA'
      meta_vals <- as.character(xbj_lab[[annot]][,1])
      cells_to_plot <- colnames(xbj_lab)[!tolower(meta_vals) %in% c("unknown", "na") & !is.na(meta_vals)]
      
      # Only plot if we have annotated cells remaining
      if (length(cells_to_plot) > 50) {
        
        # Temporary subset for the figure
        plot_obj <- subset(xbj_lab, cells = cells_to_plot)
        
        p <- DimPlot(plot_obj, 
                     reduction = "umap.mnn", 
                     group.by = annot, 
                     label = TRUE, 
                     repel = TRUE, 
                     raster = FALSE,
                     shuffle = TRUE) + # Shuffle ensures colors are mixed well
             ggtitle(paste0("Integrated (Known Only) | Lab: ", current_lab, " | ", annot)) +
             theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(paste0("Int_Clean_", current_lab, "_", annot, ".png"), plot = p, width = 12, height = 10)
        
        # Cleanup plot_obj immediately
        rm(plot_obj); gc()
        
      } else {
        print(paste0("Skipping ", annot, " for ", current_lab, ": No known annotations found."))
      }
    }
  }

  rm(xbj_lab); gc()
}