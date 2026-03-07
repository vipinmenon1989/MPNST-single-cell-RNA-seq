library(BPCells)

#' @title Open a BPCells Matrix Directory
#' @description Opens an existing BPCells matrix directory from the specified path.
#' @param path The file system path to the BPCells matrix directory.
#' @return A BPCells matrix object representing the opened directory.
mat <- open_matrix_dir("/autofs/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/counts_bpcells")

#' @title Write BPCells Matrix to HDF5
#' @description Writes the BPCells matrix object to an HDF5 file.
#'   This function uses the 'path' argument for the output file,
#'   aligning with the latest BPCells package versions.
#' @param mat The BPCells matrix object to be written.
#' @param path The output file path for the HDF5 file.
#' @param group The group name within the HDF5 file where the matrix will be stored.
#' @param compress A logical value indicating whether to compress the data in the HDF5 file.
#' @return No explicit return value; the function writes the matrix data to the specified HDF5 file.
write_matrix_hdf5(
  mat = mat, 
  path = "counts.h5", 
  group = "counts",
  compress = TRUE
)

#' @title Export Seurat Object Metadata
#' @description Exports the metadata from a Seurat object to a CSV file.
#'   This step is crucial for ensuring metadata is accessible for Python-based downstream analysis.
#' @param x The metadata data frame from the Seurat object (e.g., `seurat_obj@meta.data`).
#' @param file The name of the output CSV file, typically "metadata.csv".
#' @param row.names A logical value indicating whether to write row names (e.g., cell identifiers) to the file.
#' @return No explicit return value; the metadata is written to "metadata.csv".
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)

#' @title Export Seurat Object Gene Names
#' @description Exports the gene names (row names) from a Seurat object to a TSV file.
#'   This facilitates importing gene information into Python-based analysis environments.
#' @param x The vector of gene names (e.g., `rownames(seurat_obj)`).
#' @param file The name of the output TSV file, typically "genes.tsv".
#' @return No explicit return value; the gene names are written to "genes.tsv".
write(rownames(seurat_obj), "genes.tsv")

#' @title Export Seurat Object Barcode Names
#' @description Exports the barcode names (column names) from a Seurat object to a TSV file.
#'   This ensures cell barcode information is available for Python-based analysis.
#' @param x The vector of barcode names (e.g., `colnames(seurat_obj)`).
#' @param file The name of the output TSV file, typically "barcodes.tsv".
#' @return No explicit return value; the barcode names are written to "barcodes.tsv".
write(colnames(seurat_obj), "barcodes.tsv")

print("--- R EXPORT SUCCESSFUL: counts.h5 created ---")
