#' @title Export Single-cell RNA-seq Data to Standard Formats
#' @description This script facilitates the export of single-cell RNA-seq data from BPCells H5 format and Seurat objects
#'   into widely compatible standard formats such as Matrix Market (.mtx), CSV, and TSV files.
#' @details The script performs the following main steps:
#'   \enumerate{
#'     \item Opens a BPCells H5 matrix for raw count data.
#'     \item Streams the raw counts directly to a Matrix Market (.mtx) file.
#'     \item Loads a pre-existing Seurat object to extract and export metadata, gene names, and barcode names
#'       into separate CSV and TSV files.
#'   }
#'   This process ensures that the data is accessible for downstream analysis in various environments,
#'   particularly Python-based tools that often utilize these standard formats.
#' @importFrom BPCells open_matrix_hdf5
#' @importFrom Matrix writeMM
#' @importFrom utils write.csv
#' @importFrom base readRDS write
#' @return This script does not return any R objects. It writes several files to disk:
#'   \itemize{
#'     \item \code{counts.mtx}: A Matrix Market file containing the raw count matrix.
#'     \item \code{metadata.csv}: A CSV file containing the cell metadata from the Seurat object.
#'     \item \code{genes.tsv}: A TSV file containing the feature (gene) names.
#'     \item \code{barcodes.tsv}: A TSV file containing the cell barcode names.
#'   }
#' @note Before running this script, ensure the following files are present in the working directory:
#'   \itemize{
#'     \item \code{counts.h5}: The input BPCells H5 matrix file (with a "counts" group).
#'     \item \code{seurat_obj_pre_integration.rds}: The input Seurat object RData file.
#'   }
#' @examples
#' \dontrun{
#' # To execute this script, ensure required files are in place and then source it:
#' # source("path/to/your_export_script.R")
#' # The script will print a success message upon completion.
#' }
#' @export
library(BPCells)
library(Matrix)

# 1. Open your H5 matrix
mat <- open_matrix_hdf5(path = "counts.h5", group = "counts")

# 2. Stream to Matrix Market (.mtx)
# This writes directly to disk and is the most robust bridge to Python.
writeMM(obj = mat, file = "counts.mtx")

# 3. Export Metadata (Ensuring they are fresh)
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- EXPORT SUCCESSFUL: counts.mtx created ---")
