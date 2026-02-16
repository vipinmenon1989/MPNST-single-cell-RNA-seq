library(Seurat)
library(ggplot2)

# 1. PREPARE LAB LIST
print("--- FETCHING LAB LIST ---")
temp_obj <- readRDS("seurat_obj_pre_integration.rds")
all_labs <- unique(as.character(temp_obj$lab))
all_labs <- all_labs[!is.na(all_labs)]
rm(temp_obj); gc()

# 2. UNINTEGRATED LOOP
for (current_lab in all_labs) {
  print(paste0("--- UNINTEGRATED: ", current_lab, " ---"))
  
  # Load, Subset, and EVICT parent immediately
  pbj <- readRDS("seurat_obj_pre_integration.rds")
  pbj_lab <- subset(pbj, subset = lab == current_lab)
  rm(pbj); gc()

  # Skip labs with insufficient cells to form a UMAP
  if (ncol(pbj_lab) < 50) {
    print("Skipping: Too few cells.")
    rm(pbj_lab); gc(); next
  }

  # Processing (Local PCA/UMAP)
  pbj_lab <- NormalizeData(pbj_lab, verbose = FALSE)
  pbj_lab <- FindVariableFeatures(pbj_lab, verbose = FALSE)
  pbj_lab <- ScaleData(pbj_lab, verbose = FALSE)
  pbj_lab <- RunPCA(pbj_lab, verbose = FALSE)
  pbj_lab <- RunUMAP(pbj_lab, dims = 1:30, verbose = FALSE)

  # Plotting
  for (annot in c("wu_annot", "suppiah_annot")) {
    if (annot %in% colnames(pbj_lab@meta.data)) {
      p <- DimPlot(pbj_lab, reduction = "umap", group.by = annot, label = TRUE, repel = TRUE, raster = FALSE) +
           ggtitle(paste0("Unintegrated | Lab: ", current_lab, " | ", annot))
      ggsave(paste0("Unint_", current_lab, "_", annot, ".png"), plot = p, width = 12, height = 10)
    }
  }
  
  rm(pbj_lab); gc()
}
