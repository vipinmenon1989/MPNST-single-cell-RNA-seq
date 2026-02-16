library(Seurat)

# 1. LOAD FILE
print("--- READING RDS FILE ---")
xbj <- readRDS("seurat_obj_pre_integration.rds")

# 2. CHECK REDUCTIONS
print("--- AVAILABLE REDUCTIONS ---")
print(names(xbj@reductions))

# 3. VERDICT
if ("umap.mnn" %in% names(xbj@reductions)) {
  print("✅ SUCCESS: 'umap.mnn' IS saved in the object.")
} else if ("umap" %in% names(xbj@reductions)) {
  print("⚠️ WARNING: A UMAP exists, but it is named 'umap', not 'umap.mnn'.")
} else {
  print("❌ FAILURE: NO UMAP reduction found. You must run RunUMAP() again.")
}
