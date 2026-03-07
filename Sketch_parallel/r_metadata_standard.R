#' @title Export Seurat Object for BPCells Processing
#' @description This script loads a Seurat object, extracts its count matrix, metadata, genes, and barcodes,
#'   and exports them to disk in formats compatible with the BPCells package.
#' @details The primary count matrix from the "RNA" assay's "counts" layer is saved into a directory
#'   named "counts_bpcells" using `BPCells::write_matrix_dir`. This method optimizes storage and
#'   access for large single-cell datasets. Metadata is saved as a CSV file, while gene names
#'   and cell barcode names are saved into separate tab-separated value files.
#'   This preparation facilitates efficient downstream analysis leveraging the BPCells package's capabilities.
#' @section Output Files:
#'   - `counts_bpcells/`: A directory containing the BPCells-formatted sparse count matrix.
#'   - `metadata.csv`: A CSV file with the `seurat_obj@meta.data` table.
#'   - `genes.tsv`: A TSV file listing all gene identifiers (rownames of the Seurat object).
#'   - `barcodes.tsv`: A TSV file listing all cell barcode identifiers (colnames of the Seurat object).
#' @seealso [Seurat::readRDS()], [BPCells::write_matrix_dir()]
#' @keywords data-export Seurat BPCells single-cell
#' @author Your Name/Organization Here (Placeholder)
#' @examples
#' \dontrun{
#' # Ensure 'seurat_obj_pre_integration.rds' is present in the working directory.
#' #
#' # To run this script:
#' # source("path/to/this/script.R")
#' #
#' # After execution, the output files will be created in the current working directory.
#' # You can then inspect them:
#' # list.files()
#' # read.csv("metadata.csv")
#' }
#'
library(Seurat)
library(BPCells)

#' # 1. LOAD DATA
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")

#' # 2. WRITE MATRIX TO DISK (BPCells handles v5 layers automatically)
#' # This creates a folder 'counts_bpcells' with the 1M cell matrix
write_matrix_dir(
  mat = seurat_obj[["RNA"]]$counts,
  dir = "counts_bpcells",
  overwrite = TRUE
)

#' # 3. EXPORT METADATA & GENES
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- R EXPORT COMPLETE: Matrix stored in 'counts_bpcells' ---")
