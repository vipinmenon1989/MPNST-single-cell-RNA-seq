#' @title Check for 'umap.mnn' reduction in a Seurat object.
#' @description This script loads a Seurat object from a specified RDS file and checks for the presence of a 'umap.mnn' reduction.
#' @details It first attempts to load a Seurat object from "seurat_obj_pre_integration.rds".
#'   It then prints all available reductions in the loaded object.
#'   Finally, it provides a verdict based on the presence and naming of UMAP reductions:
#'   - "SUCCESS" if 'umap.mnn' is found.
#'   - "WARNING" if 'umap' (a different UMAP name) is found instead of 'umap.mnn'.
#'   - "FAILURE" if no UMAP reduction is found.
#' @keywords Seurat UMAP reduction check script
#' @noRd
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
