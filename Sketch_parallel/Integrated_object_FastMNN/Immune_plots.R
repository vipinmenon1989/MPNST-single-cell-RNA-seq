library(Seurat)
library(ggplot2)

# 1. GET LAB LIST AND CLEAR RAM
# -------------------------------------------------------------------------
print("--- FETCHING LAB LIST ---")
# Load briefly to get names, then delete to free ~250GB
temp_obj <- readRDS("seurat_obj_pre_integration.rds")
all_labs <- unique(as.character(temp_obj$lab))
all_labs <- all_labs[!is.na(all_labs)]
rm(temp_obj)
gc() 

# 2. SEQUENTIAL PROCESSING LOOP
# -------------------------------------------------------------------------
for (current_lab in all_labs) {
  
  print(paste0("========== STARTING LAB: ", current_lab, " =========="))
  
  # --- PART A: UNINTEGRATED (Load, Subset, Discard Parent) ---
  print("--- Loading Unintegrated Object ---")
  pbj <- readRDS("seurat_obj_pre_integration.rds")
  pbj_lab <- subset(pbj, subset = lab == current_lab)
  rm(pbj) # CRITICAL: Remove the 1M-cell parent immediately
  gc()

  print(paste0("Processing Unintegrated UMAP for: ", current_lab))
  pbj_lab <- NormalizeData(pbj_lab, verbose = FALSE)
  pbj_lab <- FindVariableFeatures(pbj_lab, verbose = FALSE)
  pbj_lab <- ScaleData(pbj_lab, verbose = FALSE)
  pbj_lab <- RunPCA(pbj_lab, verbose = FALSE)
  pbj_lab <- RunUMAP(pbj_lab, dims = 1:30, verbose = FALSE)

  # Save Unintegrated Plots
  for (annot in c("wu_annot", "suppiah_annot")) {
    if (annot %in% colnames(pbj_lab@meta.data)) {
      p1 <- DimPlot(pbj_lab, reduction = "umap", group.by = annot, label = TRUE, repel = TRUE, raster = FALSE) +
            ggtitle(paste0("Unintegrated | Lab: ", current_lab, " | ", annot))
      ggsave(paste0("Unint_", current_lab, "_", annot, ".png"), plot = p1, width = 12, height = 10)
    }
  }
  
  rm(pbj_lab)
  gc() # Ensure RAM is completely empty before loading next big file

  # --- PART B: INTEGRATED (Load, Subset, Discard Parent) ---
  print("--- Loading Integrated Object ---")
  xbj <- readRDS("seurat_obj_mnn_integrated_1M.rds")
  xbj_lab <- subset(xbj, subset = lab == current_lab)
  rm(xbj) # CRITICAL: Remove the 1M-cell parent immediately
  gc()

  # Save Integrated Plots (UMAP already exists)
  for (annot in c("wu_annot", "suppiah_annot")) {
    if (annot %in% colnames(xbj_lab@meta.data)) {
      p2 <- DimPlot(xbj_lab, reduction = "umap.mnn", group.by = annot, label = TRUE, repel = TRUE, raster = FALSE) +
            ggtitle(paste0("Integrated | Lab: ", current_lab, " | ", annot))
      ggsave(paste0("Int_", current_lab, "_", annot, ".png"), plot = p2, width = 12, height = 10)
    }
  }

  rm(xbj_lab)
  gc()
  print(paste0("--- FINISHED LAB: ", current_lab, " ---"))
}