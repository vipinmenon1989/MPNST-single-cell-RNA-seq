library(Seurat)

# 1. Load the object (if not already loaded)
 seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

print("========== OBJECT SUMMARY ==========")
print(seurat_obj)

print("\n========== REDUCTIONS AVAILABLE ==========")
# Lists all dimensionality reductions (PCA, UMAP, MNN, etc.) stored in the object
reductions <- Reductions(seurat_obj)
print(reductions)

# Detailed check of the specific reductions (dimensions computed for each)
for (red in reductions) {
  print(paste0("Reduction: ", red, " | Dimensions: ", ncol(Embeddings(seurat_obj, reduction = red))))
}

print("\n========== METADATA COLUMNS ==========")
# Lists all columns in the metadata slot
meta_cols <- colnames(seurat_obj@meta.data)
print(meta_cols)

print("\n========== RESOLUTION COLUMNS FOUND ==========")
# Specifically filters for clustering resolution columns (usually contain 'res')
res_cols <- grep("res", meta_cols, value = TRUE)
if (length(res_cols) > 0) {
  print(res_cols)
} else {
  print("WARNING: No columns containing 'res' found. Check your clustering step.")
}

print("\n========== ASSAYS AVAILABLE ==========")
# Checks if you have 'RNA', 'SCT', or 'integrated' assays
print(Assays(seurat_obj))

print("\n========== DEFAULT ASSAY ==========")
print(DefaultAssay(seurat_obj))
