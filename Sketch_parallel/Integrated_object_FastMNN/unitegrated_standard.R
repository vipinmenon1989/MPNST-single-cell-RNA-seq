#' @title Process and Plot Unintegrated Seurat Data by Lab
#' @description This script iterates through unique labs in a Seurat object, processes each lab's data independently,
#'   and generates UMAP plots based on specified annotations.
#' @details The script first loads a Seurat object to identify unique 'lab' identifiers.
#'   It then loops through each lab:
#'   1. Loads the full Seurat object again.
#'   2. Subsets the object to include only cells from the current lab.
#'   3. Skips labs with fewer than 50 cells.
#'   4. Performs standard Seurat processing steps (NormalizeData, FindVariableFeatures, ScaleData, RunPCA, RunUMAP)
#'      locally on the lab-specific subset.
#'   5. Generates UMAP plots for 'wu_annot' and 'suppiah_annot' (if available) and saves them as PNG files.
#'   Intermediate objects are removed using `rm()` and `gc()` to manage memory.
#' @import Seurat
#' @import ggplot2
#' @return PNG files are saved to the current working directory, named `Unint_[lab_name]_[annotation].png`.
NULL

library(Seurat)
library(ggplot2)

#' 1. PREPARE LAB LIST
#' This section identifies all unique 'lab' identifiers from the initial Seurat object.
print("--- FETCHING LAB LIST ---")
temp_obj <- readRDS("seurat_obj_pre_integration.rds")
all_labs <- unique(as.character(temp_obj$lab))
all_labs <- all_labs[!is.na(all_labs)]
rm(temp_obj); gc()

#' 2. UNINTEGRATED LOOP
#' This loop iterates through each unique lab, processing its data independently and generating plots.
for (current_lab in all_labs) {
  print(paste0("--- UNINTEGRATED: ", current_lab, " ---"))
  
  #' Load, Subset, and EVICT parent immediately
  #' Loads the full Seurat object and subsets it to include only cells from the current lab.
  pbj <- readRDS("seurat_obj_pre_integration.rds")
  pbj_lab <- subset(pbj, subset = lab == current_lab)
  rm(pbj); gc()

  #' Skip labs with insufficient cells to form a UMAP
  #' Checks if the current lab has at least 50 cells; otherwise, it skips processing for this lab.
  if (ncol(pbj_lab) < 50) {
    print("Skipping: Too few cells.")
    rm(pbj_lab); gc(); next
  }

  #' Processing (Local PCA/UMAP)
  #' Performs standard Seurat workflow steps on the lab-specific subset: normalization,
  #' variable feature finding, scaling, PCA, and UMAP dimensionality reduction.
  pbj_lab <- NormalizeData(pbj_lab, verbose = FALSE)
  pbj_lab <- FindVariableFeatures(pbj_lab, verbose = FALSE)
  pbj_lab <- ScaleData(pbj_lab, verbose = FALSE)
  pbj_lab <- RunPCA(pbj_lab, verbose = FALSE)
  pbj_lab <- RunUMAP(pbj_lab, dims = 1:30, verbose = FALSE)

  #' Plotting
  #' Generates and saves UMAP plots for 'wu_annot' and 'suppiah_annot' if these annotations
  #' are present in the Seurat object's metadata.
  for (annot in c("wu_annot", "suppiah_annot")) {
    if (annot %in% colnames(pbj_lab@meta.data)) {
      p <- DimPlot(pbj_lab, reduction = "umap", group.by = annot, label = TRUE, repel = TRUE, raster = FALSE) +
           ggtitle(paste0("Unintegrated | Lab: ", current_lab, " | ", annot))
      ggsave(paste0("Unint_", current_lab, "_", annot, ".png"), plot = p, width = 12, height = 10)
    }
  }
  
  rm(pbj_lab); gc()
}
