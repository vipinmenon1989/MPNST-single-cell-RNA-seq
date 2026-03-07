#' @title Explore a Seurat Object
#' @description This script loads a Seurat object from an RDS file and provides a summary of its key features,
#' including available reductions, metadata columns, clustering resolutions, and assays.
#' @details It's designed to give a quick overview of a processed Seurat object's structure
#' without performing any modifications.
#' @return Prints various summaries to the console.
#' @keywords Seurat, object exploration, single-cell RNA-seq
#' @seealso \code{\link[Seurat]{Seurat-class}}, \code{\link[Seurat]{Reductions}}, \code{\link[Seurat]{Assays}}, \code{\link[Seurat]{DefaultAssay}}
#' @export
library(Seurat)

# 1. Load the object (if not already loaded)
#' @section Load Seurat Object:
#' Loads a pre-processed Seurat object from an RDS file.
#' @param seurat_obj_clustered_1M.rds Path to the RDS file containing the Seurat object.
#' @return A Seurat object assigned to the `seurat_obj` variable.
seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

print("========== OBJECT SUMMARY ==========")
#' @section Object Summary:
#' Prints a summary of the loaded Seurat object, showing basic statistics like
#' number of cells, features, and assays.
print(seurat_obj)

print("\n========== REDUCTIONS AVAILABLE ==========")
# Lists all dimensionality reductions (PCA, UMAP, MNN, etc.) stored in the object
#' @section Available Reductions:
#' Identifies and lists all dimensionality reduction techniques (e.g., PCA, UMAP)
#' that have been computed and stored within the Seurat object.
reductions <- Reductions(seurat_obj)
print(reductions)

# Detailed check of the specific reductions (dimensions computed for each)
#' @section Reduction Dimensions:
#' Iterates through each identified reduction and reports the number of dimensions
#' computed for that specific reduction. This helps confirm the depth of each embedding.
for (red in reductions) {
  print(paste0("Reduction: ", red, " | Dimensions: ", ncol(Embeddings(seurat_obj, reduction = red))))
}

print("\n========== METADATA COLUMNS ==========")
# Lists all columns in the metadata slot
#' @section Metadata Columns:
#' Lists all column names available in the Seurat object's metadata slot.
#' These columns typically contain cell-specific annotations like cluster assignments,
#' sample origin, or QC metrics.
meta_cols <- colnames(seurat_obj@meta.data)
print(meta_cols)

print("\n========== RESOLUTION COLUMNS FOUND ==========")
# Specifically filters for clustering resolution columns (usually contain 'res')
#' @section Clustering Resolution Columns:
#' Filters the metadata columns to identify those likely containing clustering
#' resolution results (e.g., `RNA_snn_res.0.5`, `seurat_clusters`).
#' It checks for columns containing the string "res".
res_cols <- grep("res", meta_cols, value = TRUE)
if (length(res_cols) > 0) {
  print(res_cols)
} else {
  print("WARNING: No columns containing 'res' found. Check your clustering step.")
}

print("\n========== ASSAYS AVAILABLE ==========")
# Checks if you have 'RNA', 'SCT', or 'integrated' assays
#' @section Available Assays:
#' Lists all assays (e.g., "RNA", "SCT", "integrated") present in the Seurat object.
#' Assays represent different data modalities or processing pipelines.
print(Assays(seurat_obj))

print("\n========== DEFAULT ASSAY ==========")
#' @section Default Assay:
#' Reports the currently set default assay for the Seurat object. The default
#' assay is used by many Seurat functions if not explicitly specified.
print(DefaultAssay(seurat_obj))
