#' @file export_seurat_metadata.R
#' @title Export Metadata and Coordinates from a Seurat Object
#' @description This script loads a pre-integration Seurat object from an RDS file
#'   and exports its metadata, gene names, and cell barcode names to separate
#'   CSV and TSV files. This process is designed to specifically extract these
#'   components, potentially mitigating issues like 32-bit integer errors or
#'   memory constraints that can occur when handling very large Seurat objects.
#' @details The script expects the input file "seurat_obj_pre_integration.rds"
#'   to be in the current working directory. It generates "metadata.csv",
#'   "genes.tsv", and "barcodes.tsv" in the same directory.
#' @author Your Name/Organization
#' @seealso \code{\link[Seurat]{Seurat-class}}, \code{\link[base]{readRDS}}, \code{\link[utils]{write.csv}}, \code{\link[base]{write}}
#' @keywords Seurat metadata export coordinates
#' @import Seurat
#' @md
library(Seurat)

# 1. RELOAD THE OBJECT (Just for metadata)
print("--- RELOADING SEURAT OBJECT FOR METADATA ---")
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")

# 2. EXPORT METADATA & COORDINATES
# These are small files and won't trigger the 32-bit integer error
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- METADATA EXPORT COMPLETE ---")
