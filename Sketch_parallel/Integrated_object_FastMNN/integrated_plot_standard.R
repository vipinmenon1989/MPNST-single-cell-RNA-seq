#' @title Generate Clean UMAP Plots for Integrated Seurat Object by Lab
#' @description This script processes an integrated Seurat object to generate UMAP plots for each unique 'lab' group.
#' It specifically filters out cells annotated as 'unknown' or 'NA' from specified annotation columns ("wu_annot", "suppiah_annot")
#' before plotting to provide "clean" visualizations.
#' @details The script loads a pre-processed Seurat object, iterates through unique 'lab' identifiers,
#' subsets the object for each lab, and then generates UMAP plots for "wu_annot" and "suppiah_annot"
#' columns after removing 'unknown' and 'NA' entries. Plots are saved as PNG files.
#' Aggressive memory management (rm() and gc()) is used throughout the script to handle large Seurat objects.
#' @param seurat_obj_mnn_integrated_1M.rds (input file) A Seurat object containing integrated data,
#' including 'lab', 'wu_annot', and 'suppiah_annot' metadata columns, and a 'umap.mnn' reduction.
#' This file is expected to be in the current working directory.
#' @return PNG files are saved to the working directory, named `Int_Clean_[lab]_[annotation].png`.
#' Each plot displays a UMAP reduction for a specific lab and annotation, excluding 'unknown'/'NA' cells.
#' @seealso \code{\link[Seurat]{DimPlot}}, \code{\link[Seurat]{subset.Seurat}}
#' @import Seurat
#' @import ggplot2
#' @examples
#' # This script is designed to be run directly.
#' # Ensure 'seurat_obj_mnn_integrated_1M.rds' is in the working directory.
#' # To execute the script, simply source it:
#' # source("path/to/this/script.R")
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
